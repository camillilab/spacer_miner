"""
SPACER_MINER.PY
VERSION 0.5
Python 3.7+

Author: Jake Bourgeois
Email: jacob.bourgeois@tufts.edu
Affiliation: Camilli Lab, Tufts University
Date: 08-06-2020
License: BSD-3

Python script for the automated mining of spacers between CRISPR repeats and determination of protospacer/PAMs. This
program takes a BLAST XML file where the repeat sequence was BLASTed, and using those coordinates extracts spacers.
Spacers are validated by detection of complete, intact repeats on both sides. Then, each spacer is blasted against
the nr database to detect putative protospacers. Gene annotation is automatically extracted and applied to each
protospacer instance.

This script requires the Biopython library and BLASTn command line tools installed to function properly. BLASTn on
the local machine is used to help detect if detected protospacers are simply within arrays in other organisms.

If you use this script in your publication, please cite us:


Possible improvements to consider:

- PAM detection of this version is currently only implemented directly adjacent to protospacers. Might be nice
to have more control over that

- Multiple threads for NCBI web blasts to speed downstream analysis

- Allow user-customizable stringency on spacer validation

"""

from Bio import Entrez
import os
import csv
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from urllib.error import HTTPError
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline as bioblastn
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, generic_nucleotide
import time
from collections import defaultdict
import subprocess
import random
import argparse
import glob


# class container for Spacer hits
class SpacerHit:

    def __init__(self, hsp_payload, query_cov, similarity, start_coordinate_correction=0, end_coordinate_correction=0):

        self.id = hsp_payload.hit_id
        self.gid = self.id.split('|')[1]
        self.accession = self.id.split('|')[3].split('.')[0]
        self.coordinates = hsp_payload.hit_range
        self.strand = hsp_payload.hit_strand
        self.description = hsp_payload.hit_description
        self.query_coverage = query_cov
        self.query_similarity = similarity
        self.start_coordinate_correction = start_coordinate_correction  # Helps us nail the PAM if the query_cov < 1.0
        self.end_coordinate_correction = end_coordinate_correction  # Helps us nail the protospacer if query_cov < 1.0

        self.pam = None
        self.fiveprime_pam = None
        self.pam_intact = False
        self.gbk_path = None
        self.date_of_isolation = None
        self.location_of_isolation = None
        self.date_of_deposit = None
        self.annotation = None
        self.protospacer_seq = None
        self.seq_length = None

        self.in_existing_array = False

        return

    # method obtains info from corresponding gbk file to obtain PAM and annotation info
    def mine_info(self, data_dir, repeat, pam_length=6):

        # Load the seq from file
        gbk_file = os.path.join(data_dir, '{0}.gb'.format(self.accession))

        try:
            with open(gbk_file, 'r') as gbk_handle:
                record = SeqIO.read(gbk_handle, format='gb')

                core_features = record.features[0].qualifiers

                try:
                    core_feature_collection_date = core_features['collection_date'][0]
                except KeyError:
                    core_feature_collection_date = 'N/A'

                try:
                    core_feature_country = core_features['country'][0]
                except KeyError:
                    core_feature_country = 'N/A'

                deposit_date = record.annotations['date']

                # the PAM is located as the reverse complement of the spacer on the 3' end
                # You can set PAM length to anything really, if you want to try and determine it experimentally
                if self.strand == 1:
                    start_coordinate = self.coordinates[0] - self.start_coordinate_correction
                    end_coordinate = self.coordinates[1] + self.end_coordinate_correction
                    pam = record.seq[start_coordinate - pam_length:start_coordinate].reverse_complement()
                    fiveprime_pam = record.seq[end_coordinate:end_coordinate+pam_length].reverse_complement()
                    self.protospacer_seq = record.seq[start_coordinate:end_coordinate].reverse_complement()
                else:
                    start_coordinate = self.coordinates[0] - self.end_coordinate_correction
                    end_coordinate = self.coordinates[1] + self.start_coordinate_correction
                    pam = record.seq[end_coordinate:end_coordinate + pam_length]
                    fiveprime_pam = record.seq[start_coordinate - pam_length:start_coordinate]
                    self.protospacer_seq = record.seq[start_coordinate:end_coordinate]

                # get overlying annotations
                # to be fairrrr....this only gets annotations if the spacer lies perfectly within.
                overlying_annotations = [k for k in record.features
                                         if k.type == 'CDS'
                                         and k.location.start <= self.coordinates[0]
                                         and k.location.end >= self.coordinates[1]
                                         and len(k.location.parts) == 1]

                if overlying_annotations:
                    try:
                        overlying_annotations = [k.qualifiers['product'][0] for k in overlying_annotations][0]
                    except KeyError:
                        print("KeyError - trying to get note info...")
                        try:
                            overlying_annotations = [k.qualifiers['note'][0] for k in overlying_annotations][0]
                        except KeyError:
                            print("KeyError - no note information. Setting to None...")
                            overlying_annotations = None

                # add to object
                self.date_of_isolation = core_feature_collection_date
                self.location_of_isolation = core_feature_country
                self.date_of_deposit = deposit_date
                self.pam = pam
                self.fiveprime_pam = fiveprime_pam
                #if pam[0:2] == 'TT':
                #    self.pam_intact = True
                self.seq_length = len(record.seq)

                # check to see if I'm in an array
                self.in_existing_array = self._am_i_in_an_array(parent_seq=record.seq, repeat=repeat)

                # If the overlying annotations are null or contain 'hypothetical' or 'uncharacteri[cz]ed',
                # try a bit harder
                if (not overlying_annotations) or ('hypothetical' in overlying_annotations) or (
                        'uncharacter' in overlying_annotations):

                    if self.in_existing_array:
                        overlying_annotations = 'Part of CRISPR array'
                    else:

                        # print("No good annotation found...digging a bit deeper...")

                        all_cds = [k for k in record.features if k.type == 'CDS']
                        if not all_cds:
                            # there is no CDS information in the whole file!
                            overlying_annotations = 'No CDS info in file'
                        else:
                            upcheck, downcheck = None, None

                            upstream_annotations = [k for k in all_cds if k.location.end < self.coordinates[0]]
                            downstream_annotations = [k for k in all_cds if k.location.start > self.coordinates[1]]

                            cds_max_limit = 10  # Maximum number of genes to search

                            cds_search_limit_up = min(cds_max_limit, len(upstream_annotations))
                            cds_search_limit_down = min(cds_max_limit, len(downstream_annotations))
                            i = 0
                            while i < cds_search_limit_up:
                                cds = upstream_annotations[::-1][i]
                                try:
                                    up_cds = cds.qualifiers['product'][0]
                                    if ('hypothetical' not in up_cds) and ('uncharacteri' not in up_cds):
                                        # print("Found alternative upstream CDS!")
                                        upcheck = up_cds
                                        i = cds_search_limit_up
                                except KeyError:
                                    pass
                                i += 1

                            i = 0
                            while i < cds_search_limit_down:
                                cds = downstream_annotations[i]
                                try:
                                    down_cds = cds.qualifiers['product'][0]
                                    if ('hypothetical' not in down_cds) and ('uncharacteri' not in down_cds):
                                        # print("Found alternative downstream CDS!")
                                        downcheck = down_cds
                                        i = cds_search_limit_down
                                except KeyError:
                                    pass
                                i += 1

                            if upcheck or downcheck:
                                overlying_annotations = 'Target: {2}, Upstream: {0}, Downstream: {1}'.format(upcheck, downcheck, overlying_annotations)
                            else:
                                overlying_annotations = '{0} - No nearby annotated ORFs'.format(overlying_annotations)

                self.annotation = overlying_annotations
        except FileNotFoundError:
            print("Something went horribly wrong. Tried to find {0} and couldn't find it.".format(gbk_file))

        return

    # method attempts to determine if SpacerHit is within a CRISPR array or not
    # maybe check the left and right to see if the repeat exists?
    def _am_i_in_an_array(self, parent_seq, repeat):

        if self.strand == 1:
            start_coordinate = self.coordinates[0] - self.start_coordinate_correction
            end_coordinate = self.coordinates[1] + self.end_coordinate_correction
        else:
            start_coordinate = self.coordinates[0] - self.end_coordinate_correction
            end_coordinate = self.coordinates[1] + self.start_coordinate_correction

        left_seq = parent_seq[start_coordinate - len(repeat): start_coordinate]
        right_seq = parent_seq[end_coordinate: end_coordinate + len(repeat)]

        # edit 2020-08-03 - now check arrays by blasting against the repeat to catch polymorphisms and alternative array
        left_result = seq_blast(left_seq, 'repeat', q_cov=0.90, sim=0.93)
        right_result = seq_blast(right_seq, 'repeat', q_cov=0.90, sim=0.93)

        if left_result or right_result:
            # print("Repeat sequence detected: sequence {0} in hit {1} is most likely an array".format(self.accession, self.protospacer_seq))
            # print("LEFT SEQ: {0}\nRIGHT SEQ: {1}".format(left_seq, right_seq))
            return True

        return False


# class container for spacers
class Spacer:

    def __init__(self, coordinates, content, data_acc, spacer_num):

        self.coordinates = coordinates
        self.content = content
        self.spacer_accession = data_acc
        self.spacer_number = spacer_num
        self.name = '{0}_{1}'.format(self.spacer_accession, self.spacer_number)

        if self.content[0] == 'T' or self.content[0] == 'C':
            self.pyrimidine = True
        else:
            self.pyrimidine = False

        self.length = len(content)

        self.blast_hsps = []
        # self.box_spacer_blast_result = None
        return

    # method obtains information about the HSPs that BLASTed to the spacer, including PAM matching and annotation
    def download_spacer_hsps_info(self, data_dir):

        # get genbank infomation for all spacer hsps
        spacer_hit_accessions = [k.accession for k in self.blast_hsps if '{0}.gb'.format(k.accession) not in os.listdir(data_dir)]
        if spacer_hit_accessions:
            print("Taking a small breath for Entrez's servers...")
            time.sleep(60)
            get_genbank_info(accessions=spacer_hit_accessions, download_dir=data_dir, download_type='gbwithparts')
        return

    # method mines gbk files corresponding to spacer hsp hits to get info
    def mine_spacer_hsps_info(self, data_dir, repeat, pam_length):

        for hsp in self.blast_hsps:
            print("Mining for HSP information for {0}...".format(hsp.id))
            hsp.mine_info(data_dir, repeat, pam_length)
        return

    # method blasts against the known Box 2016 spacers to see if it is unique to this experiment.
    def blast_against_box_spacers(self, blast_db):

        # retrieve any perfect hit
        hit = seq_blast(d=self.content, dbname=blast_db, q_cov=0.9, sim=0.9)
        if hit:
            self.box_spacer_blast_result = hit.blast_id
        else:
            self.box_spacer_blast_result = 'None detected'

        return


# class container for CRISPR-containing bacterial hits
class CRISPRHit:

    def __init__(self, payload, format='hits'):

        # Assign attributes using the HSPFragment input
        if format == 'hits':
            self._get_hsp_info(payload)

        # Assign attributes using XML input
        if format == 'xml':
            self._get_xml_info(payload)

        # Other attributes to be defined later
        self.seq = None
        self.number_validated_repeats = 0
        self.validated_hsps = []
        self.spacers = []
        self.date_of_isolation = None
        self.location_of_isolation = None
        self.date_of_deposit = None
        self.strain = None
        self.validated_repeat_ranges = []
        return

    def _get_hsp_info(self, hit_list):

        self.description = hit_list.description
        self.accession = hit_list.accession
        self.hsps = [k for k in hit_list.hsps]
        self.number_hits = len(self.hsps)
        self.strand = self.hsps[0].hit_strand

        return

    def _get_xml_info(self, xml):
        # do stuff
        return

    # method finds spacers for class hsps
    def find_spacers(self):

        # using hsps ranges, obtain spacers - ensure we are on correct strand!
        # find spacer coordinates
        hit_ranges = []
        for hit in self.hsps:
            hit_start = hit.hit_start
            hit_end = hit.hit_end
            hit_ranges.append((hit_start, hit_end))

        hit_ranges = sorted(hit_ranges)

        for i in range(1, len(hit_ranges)-2):
            current_hit_range = hit_ranges[i]
            next_hit_range = hit_ranges[i+1]

            # Looking at the original array in O395, the spacer repeat preceeds the spacer information at the start of
            # the array, and follows the final spacer in the array. To grab perfect spacers, we need to ensure
            # that we only report spacers flanked by perfect repeats.


            """
            if i == 0 and self.strand == -1:
                try:
                    spacer_coordinate = (current_hit_range[0]-33, current_hit_range[0])
                    spacer_content = self.seq[spacer_coordinate[0]:spacer_coordinate[1]].reverse_complement()
                    self.spacers.append(Spacer(coordinates=spacer_coordinate, content=spacer_content, data_acc=self.accession))
                except IndexError:
                    print('IndexError!')
            """

            try:
                spacer_coordinate = (current_hit_range[1], next_hit_range[0])
                spacer_content = self.seq[spacer_coordinate[0]:spacer_coordinate[1]]
                spacer_num = len(hit_ranges) - 2 - i

                # If we are on the reverse strand, then the order of the spacers is reversed
                # Here, the most historic spacer is deemed spacer #1, and so forth
                if self.strand == -1:
                    spacer_content = spacer_content.reverse_complement()
                    spacer_num = i

                # If there's some gunk in the spacer or it's absurdly long, toss it.
                if 'NN' in spacer_content or len(spacer_content) > 200:
                    print("Abnormal spacer detected! {0}".format(spacer_coordinate))

                # If the repeats aren't in the validated repeat set, ignore
                elif (current_hit_range not in self.validated_repeat_ranges) or (next_hit_range not in self.validated_repeat_ranges):
                    print("Spacer lies between imperfect repeats!")
                    print("Left Repeat: {0}".format(self.seq[current_hit_range[0]:current_hit_range[1]]))
                    print("Right Repeat: {0}".format(self.seq[next_hit_range[0]:next_hit_range[1]]))
                    # time.sleep(5)

                else:

                    self.spacers.append(Spacer(coordinates=spacer_coordinate, content=spacer_content, data_acc=self.accession,
                                               spacer_num=spacer_num))

            except IndexError:
                print("Unable to get spacer coordinate - sequence index error")

        # finally, to make the data look pretty, reverse the internal order of spacers if we are on the positive strand
        if self.strand == 1:
            self.spacers = self.spacers[::-1]

        return

    # method validates repeats using given criteria, such as perfect matching
    def validate_repeats(self, repeat_seq):
        for k in self.hsps:

            # get repeat seq and query seq for each hsp
            hit_seq = k.hit.seq

            # adhere to criteria
            if repeat_seq == hit_seq:
                self.validated_hsps.append(k)
                hit_start = k.hit_start
                hit_end = k.hit_end
                self.validated_repeat_ranges.append((hit_start, hit_end))

            # overwrite hsps and number_hits, or save to a new validated set?
            self.number_validated_repeats = len(self.validated_hsps)

        return

    # method blasts spacers to see if there are novel targets
    def blast_spacers(self, time_between_queries=10, time_between_attempts=300, max_attempts=3, save_dir=os.getcwd(), database='nr'):

        # blast spacers against refseq? viral genomes?

        # Create FASTA record of all spacers
        test_records = [SeqRecord(seq=k.content, id=k.name) for k in self.spacers]

        if len(test_records) < 1:
            print("Error! Accession {0} has no valid spacers. Skipping...".format(self.accession))
        else:

            with open('spacers_record.fa', 'w') as test_handle:
                SeqIO.write(test_records, test_handle, format='fasta')

            # Load record into memory
            with open('spacers_record.fa', 'r') as test_handle:
                query = test_handle.read()

            print("Waiting {0} seconds for next query...".format(time_between_queries))
            time.sleep(time_between_queries)
            print('done.')
            blast_record_dir = os.path.join(save_dir, '{0}.xml'.format(self.accession))
            attempts = 0
            while attempts < max_attempts:
                print("BLASTing {0} to NR (attempt {1} of {2})....".format(self.accession, attempts, max_attempts))
                try:
                    with NCBIWWW.qblast(program='blastn', database=database, sequence=query) as blast_handle:
                        blast_result = blast_handle.read()
                        with open(blast_record_dir, 'w') as output_handle:
                            output_handle.write(blast_result)
                        print("done. Written to {0}".format(blast_record_dir))
                        attempts = max_attempts
                except HTTPError as e:
                    print("\nHTTP ERROR! {0}".format(e))
                    attempts += 1
                    print("Waiting {0} seconds...".format(time_between_attempts))
                    time.sleep(time_between_attempts)
                except ValueError as e:
                    print("\nVALUE ERROR! {0}".format(e))
                    attempts += 1
                    print("Waiting {0} seconds...".format(time_between_attempts))
                    time.sleep(time_between_attempts)

        return

    # method appends sequence data to hit
    def add_seq(self, payload):
        with open(payload, 'r') as input_handle:
            goods = SeqIO.read(input_handle, format='fasta')
            self.seq = goods.seq
        return

    # method adds genbank info to hit
    def add_genbank_info(self, payload):
        with open(payload, 'r') as input_handle:
            goods = SeqIO.read(input_handle, format='gb')

            core_features = goods.features[0].qualifiers

            try:
                core_feature_collection_date = core_features['collection_date'][0]
            except KeyError:
                core_feature_collection_date = 'N/A'

            try:
                core_feature_country = core_features['country'][0]
            except KeyError:
                core_feature_country = 'N/A'

            deposit_date = goods.annotations['date']

            self.seq = goods.seq
            self.date_of_isolation = core_feature_collection_date
            self.location_of_isolation = core_feature_country
            self.date_of_deposit = deposit_date
            self.description = goods.description

    # method obtains spacer information from a accession BLAST XML query
    def mine_spacer_info(self, blast_results_dir, query_cov, similarity):

        # Load the proper blast XML
        blast_xml = os.path.join(blast_results_dir, '{0}.xml'.format(self.accession))

        # Parse the XML
        print("Parsing {0}...".format(blast_xml))
        with open(blast_xml, 'r') as xml_handle:
            for blast_result in SearchIO.parse(xml_handle, format='blast-xml'):

                # Obtain hit information
                spacer_hits = blast_result.hsps
                spacer_id = blast_result.id
                print("Obtained {0} hits.".format(len(spacer_hits)))

                # Find the spacer within this class using a quick n dirty list comprehension
                try:
                    bug_spacer = [k for k in self.spacers if k.name == spacer_id][0]

                    # Get hits that fit the profile
                    i = 0
                    for hit in spacer_hits:
                        # hit_query_coverage = (hit.query_range[1] - hit.query_range[0]) / blast_result.seq_len
                        hit_query_coverage = hit.query_span / blast_result.seq_len
                        # hit_similarity = hit.aln_annotation['similarity'].count('|') / blast_result.seq_len
                        hit_similarity = hit.aln_annotation['similarity'].count('|') / hit.aln_span

                        if hit_query_coverage >= query_cov:
                            if hit_similarity >= similarity:

                                # This is an interesting hit, add it to the list
                                i += 1
                                # If we're not matching on the end of the blast hit, then we need to be sure to capture
                                # the entire protospacer
                                end_correction = blast_result.seq_len - hit.query_end
                                bug_spacer.blast_hsps.append(SpacerHit(hsp_payload=hit, query_cov=hit_query_coverage,
                                                                       similarity=hit_similarity,
                                                                       start_coordinate_correction=hit.query_start,
                                                                       end_coordinate_correction=end_correction))
                    print("Added {0} HSPS to spacer {1}".format(i, bug_spacer.name))

                except IndexError:
                    pass  # This error occurs when using old BLAST data for spacers. The numbering should be preserved, however
        return

    # method obtains gbk information for spacer HSP hits
    def obtain_spacer_hit_info(self, entrez_dir, repeat, pam_length):

        for spacer in self.spacers:

            # download spacer hit hsp genbank info
            print("Downloading Entrez data for spacer {0} to {1}...".format(spacer.name, entrez_dir))
            spacer.download_spacer_hsps_info(data_dir=entrez_dir)

            # mine for information
            print("Obtaining HSP information for spacer {0}...".format(spacer.name))
            spacer.mine_spacer_hsps_info(data_dir=entrez_dir, repeat=repeat, pam_length=pam_length)

            # also determine if the spacer already exists in the original Box2016 paper
            # spacer.blast_against_box_spacers(blast_db='box_spacers')

        return


# ------NON-CLASS FUNCTIONS-------
# Blasts a sequence over a local database.
def blast(seq, dbname, out_xml, prog):

    # Run the blast command and save the xml
    if prog == 'blastn':
        blast_prog = bioblastn(query=seq, db=dbname, outfmt=5, out=out_xml, num_threads=4, task='blastn-short')
        _, _ = blast_prog()
    else:
        print("Program {0} not recognized.".format(prog))
    return


# Create a local blast database of provided FASTA.
def make_local_blast_db(fasta, dbname, dbtype='nucl'):
    # dbname_path = os.path.join('/usr/local/ncbi/blast/db', dbname)
    subprocess.run(['makeblastdb', '-in', fasta, '-parse_seqids', '-title', dbname, '-dbtype', dbtype, '-out', dbname],
                   stdout=subprocess.DEVNULL)
    return


# returns highest blast hit
def seq_blast(d, dbname, q_cov=1.0, sim=1.0):

    my_record = SeqRecord(seq=d, id='MY_SEQ')
    with open('temp.fasta', 'w') as o:
        SeqIO.write(my_record, o, 'fasta')
    blast('temp.fasta', dbname=dbname, out_xml='temp.xml', prog='blastn')
    blast_result = SearchIO.read('temp.xml', 'blast-xml')
    best_hit = 'None detected'
    bh = None
    hit_query_coverage = 0
    hit_similarity = 0
    if blast_result.hits:
        bh = blast_result.hits[0]
        hsp = bh.hsps[0]
        best_hit_name = bh.blast_id
        query_range = hsp.query_range
        hit_range = hsp.hit_range

        hit_query_coverage = hsp.query_span / bh.seq_len
        hit_similarity = hsp.aln_annotation['similarity'].count('|') / hsp.aln_span

        best_hit = '{0} at {1} over query range {2}'.format(best_hit_name, hit_range, query_range)

    # print("Best hit: {0}".format(best_hit))
    if hit_query_coverage >= q_cov and hit_similarity >= sim:
        # print("Spacer detected: {0}".format(best_hit))
        return bh
    else:
        return None


# gets accession info using Entrez
def get_genbank_info(accessions, download_dir, batch_size=1000, download_type='fasta', max_attempts=3, time_between_attempts=300):

    epost_handle = Entrez.epost(db='nuccore', id=','.join(accessions))
    epost_result = Entrez.read(epost_handle)

    # obtain WebEnv and Query_Key for EFetch query
    query_key = epost_result['QueryKey']
    webenv = epost_result['WebEnv']

    # define output type
    if download_type == 'fasta':
        download_suffix = '.fa'
        write_format = 'fasta'
    elif download_type == 'gbwithparts':
        download_suffix = '.gb'
        write_format = 'gb'
    else:
        download_suffix = '.gb'
        write_format = 'gb'

    # Using efetch, grab the sequence and download it to local storage. Do this in a batch of 100

    attempts = 1
    while attempts < max_attempts:
        try:
            with Entrez.efetch(db='nuccore', retmax=batch_size, rettype=download_type, retmode='text', webenv=webenv, query_key=query_key) as efetch_handle:
                efetch_result = SeqIO.parse(efetch_handle, format=write_format)
                for record in efetch_result:
                    output_filename = os.path.join(download_dir, '{0}{1}'.format(record.name.split('.')[0], download_suffix))
                    with open(output_filename, 'w') as output_file:
                        SeqIO.write(record, output_file, format=write_format)
                        print("Wrote {0}".format(output_filename))
            attempts = max_attempts
        except HTTPError as e:
            print("\nHTTP ERROR! {0}".format(e))
            attempts += 1
            print("Waiting {0} seconds...".format(time_between_attempts))
            time.sleep(time_between_attempts)
        except ValueError as e:
            print("\nVALUE ERROR! {0}".format(e))
            attempts += 1
            print("Waiting {0} seconds...".format(time_between_attempts))
            time.sleep(time_between_attempts)

    return


# dump file 1 - overview of crisprhit from original blast information
def dump_overview_file(bugs, dumpfile):

    headers = ['Accession',
               'Description',
               'Date of Deposit',
               'Date of Isolation',
               'Location of Isolation',
               'Number of Repeats Detected',
               'Number of Validated Repeats',
               'Number of Mined Spacers']

    with open(dumpfile, 'w') as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=headers, delimiter='\t')
        writer.writeheader()

        data_dict = dict()

        # generate data field
        for bug in bugs:
            data_fields = [bug.accession,
                           bug.description,
                           bug.date_of_deposit,
                           bug.date_of_isolation,
                           bug.location_of_isolation,
                           bug.number_hits,
                           bug.number_validated_repeats,
                           len(bug.spacers)]
            i = 0
            for i in range(0, len(headers)):
                data_dict[headers[i]] = data_fields[i]

            # write payload
            writer.writerow(data_dict)

    return


# dump file 2 - overview of mined spacer information from CRISPRHits
def dump_spacer_info_file(bugs, dumpfile):

    headers = ['CRISPRHit Accession',
               'CRISPRHit Description',
               'Spacer Number',
               'Spacer Sequence',
               'Spacer Coordinate Start',
               'Spacer Coordinate End',
               'Spacer Length',
               'Number HSPS to Spacer Content']

    with open(dumpfile, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerow(headers)
        # generate data string
        for bug in bugs:
            bug_info = [bug.accession,
                        bug.description]
            for spacer in bug.spacers:
                spacer_info = [spacer.spacer_number,
                               spacer.content,
                               spacer.coordinates[0],
                               spacer.coordinates[1],
                               spacer.length,
                               len(spacer.blast_hsps)]
                data_string = bug_info + spacer_info
                writer.writerow(data_string)
    return


# dump file 3 - the whole shebang
def dump_full_file(bugs, dumpfile):

    # 2020-08-03 Added spacer sequence, protospacer reverse complement as columns

    headers = ['CRISPRHit Accession',
               'CRISPRHit Description',
               'Spacer Number',
               'Spacer Sequence',
               'Total Spacers in Array',
               'Spacer Hit Accession',
               'Spacer Hit Description',
               'Spacer Hit Annotations',
               'Spacer Hit Sequence Length',
               'Spacer Hit Date of Deposit',
               'Spacer Hit Date of Collection',
               'Spacer Hit Location of Isolation',
               'Spacer Hit Protospacer',
               'Spacer Hit Protospacer Reverse Complement',
               'Spacer Hit Start Position',
               'Spacer Hit End Position',
               'Spacer Hit Query Coverage',
               'Spacer Hit Query Similarity',
               'Spacer Hit 5\' PAM',
               'Spacer Hit 3\' PAM',
               'In existing array?']

    with open(dumpfile, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerow(headers)

        # generate data string
        for bug in bugs:
            bug_info = [bug.accession, bug.description]
            for spacer in bug.spacers:
                spacer_info = [spacer.spacer_number, spacer.content, len(bug.spacers)]
                for hsp in spacer.blast_hsps:
                    hsp_info = [hsp.accession,
                                hsp.description,
                                hsp.annotation,
                                hsp.seq_length,
                                hsp.date_of_deposit,
                                hsp.date_of_isolation,
                                hsp.location_of_isolation,
                                hsp.protospacer_seq,
                                hsp.protospacer_seq.reverse_complement(),
                                hsp.coordinates[0],
                                hsp.coordinates[1],
                                hsp.query_coverage,
                                hsp.query_similarity,
                                hsp.fiveprime_pam,
                                hsp.pam,
                                hsp.in_existing_array]
                    data_string = bug_info + spacer_info + hsp_info
                    writer.writerow(data_string)

    return


# dump file 4 - spacer histogram
def dump_spacer_histogram_file(bugs, dumpfile):

    spacer_histogram = defaultdict(int)
    spacer_parent_dict = defaultdict(set)
    spacer_match_dict = defaultdict(set)
    spacer_accession_dict = defaultdict(list)
    spacer_protospacer_dict = defaultdict(bool)

    # iterate through spacers
    for bug in bugs:
        for spacer in bug.spacers:
            content = str(spacer.content)
            spacer_histogram[content] += 1
            spacer_parent_dict[content].add(bug.accession)
            # spacer_match_dict[content].add(spacer.box_spacer_blast_result)
            spacer_accession_dict[content].append(spacer.name)

            # This little loop should only trigger if not seen before, and trip to True if
            # it encounters any matching hsp to a non-array spacer.
            if not spacer_protospacer_dict[content]:
                spacer_protospacer_dict[content] = False
                if spacer.blast_hsps:
                    for hsps in spacer.blast_hsps:
                        if not hsps.in_existing_array:
                            spacer_protospacer_dict[content] = True

    # dump histogram
    with open(dumpfile, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerow(['Spacer', 'Number Occurrences',
                         'Spacer Accessions', 'Novel Spacer?',
                         'Did Spacer match to Protospacer?'])

        for entry in sorted(spacer_histogram.items(), key=lambda x: x[1], reverse=True):
            writer.writerow([entry[0],
                             entry[1],
                             ', '.join(spacer_accession_dict[entry[0]]),
                             spacer_match_dict[entry[0]],
                             spacer_protospacer_dict[entry[0]]])
    return


# dump file 5 - file showing matched SpacerHit accessions to new arrays - could be novel source of spacers?
def dump_spacer_hit_array_file(bugs, dumpfile):

    # collect unique list (set) of elements in SpacerHit where spacers were found in arrays
    new_spacerhit_accessions = set()
    spacerhit_info = dict()

    for bug in bugs:
        for spacer in bug.spacers:
            for hsp in spacer.blast_hsps:
                if hsp.in_existing_array:
                    new_spacerhit_accessions.add(hsp.accession)
                    spacerhit_info[hsp.accession] = {'Description': hsp.description}

    with open(dumpfile, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerow(['Accession', 'Description'])

        for entry in sorted(list(new_spacerhit_accessions)):
            payload = [entry, spacerhit_info[entry]['Description']]
            writer.writerow(payload)

    return


# dump file 6 - equivalent to file 3, but demands perfect PAM and perfect query similarity - ie perfect spacers
def dump_full_perfect_file(bugs, dumpfile):

    headers = ['CRISPRHit Accession',
               'CRISPRHit Description',
               'Spacer Number',
               'Total Spacers in Array',
               'Spacer Hit Accession',
               'Spacer Hit Description',
               'Spacer Hit Annotations',
               'Spacer Hit Sequence Length',
               'Spacer Hit Date of Deposit',
               'Spacer Hit Date of Collection',
               'Spacer Hit Location of Isolation',
               'Spacer Hit Protospacer',
               'Spacer Hit Start',
               'Spacer Hit End',
               'Spacer Hit Query Coverage',
               'Spacer Hit Query Similarity',
               'Spacer Hit PAM',
               'In existing array?']

    with open(dumpfile, 'w') as output_handle:
        writer = csv.writer(output_handle, delimiter='\t')
        writer.writerow(headers)

        # generate data string
        for bug in bugs:
            bug_info = [bug.accession, bug.description]
            for spacer in bug.spacers:
                spacer_info = [spacer.spacer_number, len(bug.spacers)]
                for hsp in spacer.blast_hsps:
                    hsp_info = None
                    # You may choose to keep pam_intact as a condition if you so please.
                    if hsp.query_similarity == 1 and hsp.query_coverage == 1:
                        hsp_info = [hsp.accession,
                                    hsp.description,
                                    hsp.annotation,
                                    hsp.seq_length,
                                    hsp.date_of_deposit,
                                    hsp.date_of_isolation,
                                    hsp.location_of_isolation,
                                    hsp.protospacer_seq,
                                    hsp.coordinates[0],
                                    hsp.coordinates[1],
                                    hsp.query_coverage,
                                    hsp.query_similarity,
                                    hsp.pam,
                                    hsp.in_existing_array]
                    if hsp_info:
                        data_string = bug_info + spacer_info + hsp_info
                        writer.writerow(data_string)

    return


# dump file seven - potentially useful run statistics
def dump_stats_file(bugs, dumpfile):
    return

# ---------------------------MAIN-------------------------
#


def main(blast_file, repeat_sequence, cred_file, pam_length=6, qcov=0.96, qsim=0.93, sampling=False, number_accessions=0, out_prfx=''):

    blast_data_file_list = [blast_file]
    crispr_hits = []
    repeat_seq = repeat_sequence

    # number of accessions to sample
    number_of_accessions_to_sample = number_accessions

    # blast cutoffs
    query_coverage_cutoff = qcov  # 0.96 allows for one base from the query (~32-33bp) to be missing (capture 5'nt mismatches)
    query_similarity_cutoff = qsim  # 0.93 allows for up to two mismatches in the protospacer

    print("Query Coverage Cutoff: {0}\nQuery Similarity Cutoff: {1}".format(query_coverage_cutoff, query_similarity_cutoff))
    all_crispr_hit_accessions = []

    # my ncbi info - please don't share haha
    entrez_options = dict()
    with open(cred_file, 'r') as cred_handle:
        for line in cred_handle:
            entrez_options[line.split('=')[0]] = line.split('=')[1][:-1]
    if not entrez_options:
        print("No credentials detected in {0}! Please be aware this makes you susceptible to NCBI timeout and IP ban".format(cred_file))
    else:
        print('Entrez options: {0}'.format(' '.join([item for k in entrez_options for item in (k, entrez_options[k])])))
    Entrez.email = entrez_options['email']
    Entrez.api_key = entrez_options['api_key']

    for blast_data in blast_data_file_list:
        print("Using repeat sequence {0} and BLAST XML {1}".format(repeat_seq, blast_data))
        with open(blast_data, 'r') as blast_xml:
            for qresult in SearchIO.parse(blast_xml, format='blast-xml'):
                for hit in qresult.hits:
                    hit_accession = hit.accession
                    if hit_accession not in all_crispr_hit_accessions:
                        crispr_hits.append(CRISPRHit(payload=hit, format='hits'))
                        all_crispr_hit_accessions.append(hit_accession)

    # validate hits using 28bp consensus repeat, demand perfect matches
    for k in crispr_hits:
        k.validate_repeats(repeat_seq=repeat_seq)

    # prune the list of crispr_hits
    valid_crispr_hits = [k for k in crispr_hits if k.number_validated_repeats > 0]
    print("Obtained {0} valid accessions with perfect repeats from {1} BLAST results.".format(len(valid_crispr_hits), len(crispr_hits)))

    # sample
    if sampling:
        print("Sampling {0} accessions...".format(number_of_accessions_to_sample))
        valid_crispr_hits = random.sample(valid_crispr_hits, number_of_accessions_to_sample)

    print("Processing {0} accessions...".format(len(valid_crispr_hits)))

    # get the accessions
    valid_accessions = [k.accession for k in valid_crispr_hits]

    # get the genbank info

    # Directory containing Entrez Data
    data_dir = os.path.join(os.getcwd(), 'entrez_data')
    need_data_accessions = [k for k in valid_accessions if '{0}.gb'.format(k) not in os.listdir(data_dir)]
    if need_data_accessions:
        get_genbank_info(need_data_accessions, data_dir, download_type='gbwithparts')

    # Mine the spacer information from CRISPRHit files
    for bug in valid_crispr_hits:
        bug.add_genbank_info(payload=os.path.join(data_dir, '{0}.gb'.format(bug.accession)))
        bug.find_spacers()

    # Prune the list again to remove CRISPRHits that had no spacer content
    valid_spacer_hits = [k for k in valid_crispr_hits if len(k.spacers) > 0]
    print("Obtained {0} accessions with spacers from {1} validated CRISPR BLAST hits.".format(len(valid_spacer_hits), len(valid_crispr_hits)))

    # How many spacers did I find?
    i = 0
    for bug in valid_spacer_hits:
        for spacer in bug.spacers:
            i += 1
    print('{0} total spacers mined.'.format(i))

    # Directory for BLAST results
    blast_dir = os.path.join(os.getcwd(), 'blast')
    # Now let's set up the blast loop

    # DEBUG!!
    # valid_spacer_hits = valid_spacer_hits[0:1]

    # Generate local blastdb for repeat sequence
    with open('repeat.fasta', 'w') as output_handle:
        record = SeqRecord(seq=Seq(repeat_seq, alphabet=generic_nucleotide), id='repeat_seq')
        SeqIO.write(record, output_handle, format='fasta')
    make_local_blast_db(fasta='repeat.fasta', dbname='repeat')

    for bug in valid_spacer_hits:

        # Get list of spacers and BLAST to NR over the web if needed
        if '{0}.xml'.format(bug.accession) not in os.listdir(blast_dir):
            print("No BLAST data detected for {0}. BLASTing...".format(bug.accession))
            bug.blast_spacers(save_dir=blast_dir, time_between_queries=60)

        # Now get the valid hits out of the BLAST results
        print("Obtaining spacer information for {0}...".format(bug.accession))
        bug.mine_spacer_info(blast_results_dir=blast_dir, query_cov=query_coverage_cutoff, similarity=query_similarity_cutoff)

        # Obtain the genbank information for the spacer hsps
        print("Obtaining spacer HSP information for {0}...".format(bug.accession))
        bug.obtain_spacer_hit_info(entrez_dir=data_dir, repeat=repeat_seq, pam_length=pam_length)

    # Now finally, let's try to dump some data to file
    print("Writing analysis files...")

    overview_file = '{0}_crispr_hit_overview.tsv'.format(out_prfx)
    spacer_file = '{0}_spacer_overview.tsv'.format(out_prfx)
    full_file = '{0}_crispr_repeat_mined_data.tsv'.format(out_prfx)
    histo_file = '{0}_spacer_histogram.tsv'.format(out_prfx)
    array_file = '{0}_potential_new_crispr_arrays.tsv'.format(out_prfx)
    perfect_file = '{0}_perfect_crispr_repeat_mined_data.tsv'.format(out_prfx)

    print("Writing overview file to {0}...".format(overview_file))
    dump_overview_file(valid_spacer_hits, dumpfile=overview_file)

    print("Writing spacer information file to {0}...".format(spacer_file))
    dump_spacer_info_file(valid_spacer_hits, dumpfile=spacer_file)

    print("Writing protospacer information file to {0}...".format(full_file))
    dump_full_file(valid_spacer_hits, dumpfile=full_file)

    print("Writing spacer frequency file to {0}...".format(histo_file))
    dump_spacer_histogram_file(valid_spacer_hits, dumpfile=histo_file)

    print("Writing potential new CRISPR array accessions to {0}...".format(array_file))
    dump_spacer_hit_array_file(valid_crispr_hits, dumpfile=array_file)

    print("Writing curated protospacers to {0}...".format(perfect_file))
    dump_full_perfect_file(valid_spacer_hits, dumpfile=perfect_file)

    # Clean up BLASTdb files
    print("Cleaning temporary files...")
    prefixes = ['repeat']
    for prefix in prefixes:
        to_remove = glob.glob(os.path.join(os.getcwd(), '{0}*'.format(prefix)))
        for item in to_remove:
            print("Deleting {0}...".format(item))
            os.remove(item)

    return 1


if __name__ == '__main__':

    # Entry point

    # VARS
    parser = argparse.ArgumentParser()

    parser.add_argument('BLAST_XML', help='BLAST XML file of the repeat sequence coordinates')
    parser.add_argument('REPEAT_SEQUENCE', help='Repeat unit of CRISPR array')
    parser.add_argument('--query_coverage_threshold', '-c', help='Query coverage cutoff for protospacers (default=0.9)', type=float, default=0.9)
    parser.add_argument('--query_similarity_threshold', '-s', help='Query similarity cutoff for protospacers (default=0.9)', type=float, default=0.9)
    parser.add_argument('--pam_length', '-p', help='Number of nucleotides 3\' to obtain relative to the protospacer (default=6)', type=int, default=6)
    parser.add_argument('--output_prefix', '-o', help='Output file prefix.', default='', type=str)

    credentials_file = os.path.join(os.getcwd(), 'spacer_miner_creds.txt')

    if not os.path.exists(credentials_file):
        print("Credentials file {0} not found! Please place in working directory".format(credentials_file))

    args = parser.parse_args()

    # execute
    code = main(blast_file=args.BLAST_XML, repeat_sequence=args.REPEAT_SEQUENCE, pam_length=args.pam_length, qcov=args.query_coverage_threshold, qsim=args.query_similarity_threshold, cred_file=credentials_file, out_prfx=args.output_prefix)

    if code == 1:
        print("Program successfully completed!")
    else:
        print("Program completed with errors.")





