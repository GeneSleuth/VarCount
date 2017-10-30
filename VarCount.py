import sys
from Patients import Patients
from Enumerations import Genotypes, MutationEffects, Ancestry, SearchLevel
from Mutations import Mutation, Header, CodingGeneMutation, CodingGeneMutationHeader
import zipfile
import os
import time
import argparse


class CheckExome:
    """
    This class is designed to filter out intron variants.
    If your vcf file contains whole-genome variants,
    but you only want exon variants,
    and you have a .bed file identifying the location of the variants
    Pass in the chromosome number and position of a variant into check()
    """
    def __init__(self, bed_file_path):
        self.bed_file_path = bed_file_path
        self.bed_file = open(self.bed_file_path, mode='r')
        self.chromosome_number = -1
        self.chromosome_begin = -1
        self.chromosome_end = -1

    def check(self, check_chromosome_number, check_chromosome_position):
        if check_chromosome_number < self.chromosome_number or \
                (check_chromosome_position < self.chromosome_begin and check_chromosome_number == self.chromosome_number):
            self.bed_file.close()  # If you ask it to check a previous line, it has to go back and reopen
            self.bed_file = open(self.bed_file_path, mode='r')
        while self.chromosome_number != check_chromosome_number or self.chromosome_end <= check_chromosome_position:
            line = self.bed_file.readline().strip()
            self.chromosome_number, self.chromosome_begin, self.chromosome_end = map(int, line.split("\t")[0:3])
        if self.chromosome_begin < check_chromosome_position:
            return True
        else:
            return False

    def check_mutation(self, mutation):
        check_chromosome_number = int(mutation["#CHROM"])
        check_chromosome_position = int(mutation["POS"])
        return self.check(check_chromosome_number, check_chromosome_position)

    def __del__(self):
        self.bed_file.close()


class Data:
    """
    This class is a wrapper for a whole vcf file
    It reads in the input file line by line;
        It keeps track of the Header, and one Mutation at a time
    It also keeps a dictionary of coding_genes that have been affected by any of the mutations
        and constantly updates this dictionary as it sees new mutations
    It acts as a pipe for mutations, collecting their effects on coding_genes
    """
    def __init__(self, annotated_file_path, output_file_path, count_file_path, subject_info_file_path=None, bed_file_path=None,
                 valid_genotypes=None, valid_mutation_effects=None, valid_ancestries=None, search_level=None, MAF=0.02):
        self.abbreviated_titles = {}
        """
        Abbreviated titles is a dictionary that is shallow copied into all the instances of Mutation
        It's role is to provide an easy lookup for abbreviated titles that is shared between mutations
        The first time you look up mutation['Ancestral'] instead of mutation["Ancestral Allele"], the program will
            have to look through all the keys to find if any match
        Afterwards, however, an entry will be added to abbreviated_titles;
            abbreviated_titles['Ancestral'] = 'Ancestral Allele'
        Henceforth, the program will use this lookup if possible, instead of having to iterate through all the keys
        """
        self.annotated_file_path = annotated_file_path
        self.output_file_path = output_file_path
        self.subject_info_file_path = subject_info_file_path
        self.bed_file_path = bed_file_path
        self.MAF = MAF

        if self.bed_file_path is None:
            self.check_exome = None
        else:
            self.check_exome = CheckExome(self.bed_file_path)
        if valid_genotypes is None:
            self.valid_genotypes = [Genotypes.CompoundHeterozygotes]
        else:
            self.valid_genotypes = valid_genotypes
        if valid_mutation_effects is None:
            self.valid_mutation_effects = MutationEffects.default()
        else:
            self.valid_mutation_effects = valid_mutation_effects
        if valid_ancestries is None:
            self.valid_ancestries = Ancestry.all()
        else:
            self.valid_ancestries = valid_ancestries
        self.search_level = search_level

        self.mutation = None
        self.mutations = []
        self.header = None
        """
        Mutation class is a wrapper around a dictionary, and Header is a wrapper around a list.
        Header has the list of titles, in order. Mutation's dictionary is intrinsically unordered.
        To write the first column; The title is the first entry in Header's list
        The column values can be looked up in Mutation's dictionary
        """

        self.ch_counts = {}
        self.ch_count_header = None
        self.seen_ids = set()
        self.output_file = None
        self.mutation_effect_lookup = MutationEffects.str_lookup()
        with open(self.output_file_path, 'w') as self.output_file:
            if zipfile.is_zipfile(self.annotated_file_path):
                with zipfile.ZipFile(self.annotated_file_path, mode='r') as zf:
                    with zf.open(os.path.basename(self.annotated_file_path)[:-4], mode='r') as f:
                        self.parse_file_object(f)
            else:
                with open(self.annotated_file_path, 'r') as f:
                    self.parse_file_object(f)
        for coding_gene, ch_count in self.ch_counts.items():
            for patient in ch_count.header.patient_columns:
                if ch_count[patient].value is not False:
                    ch_count["COUNT"] += 1
                    for key, value in self.patients.patients[patient].data.items():
                        ch_count[key + "=" + value] += 1
        self.save_count_file(count_file_path)

    def parse_file_object(self, file_object):
        """
        Main loop of program
        """
        first_line = decode(file_object.readline())
        while first_line[:2] == "##":  # Gets rid of metadata header
            first_line = decode(file_object.readline())
        self.header = Header(first_line.strip())
        if self.subject_info_file_path is not None:
            self.patients = Patients(subject_info_file_path=self.subject_info_file_path)
        else:
            self.patients = Patients(patients=self.header.patient_columns)
        self.ch_count_header = CodingGeneMutationHeader(self.header, self.patients)
        self.output_file.write(repr(self.header))
        start = time.time()
        count = 0
        for line in file_object:
            if line == "":
                continue
            line = decode(line)
            count += 1
            if count % 1000 == 0:
                print("Line number: %s\nTime: %s" % (count, time.time() - start))
                start = time.time()
            self.mutation = Mutation(self.header, line.strip(), self.abbreviated_titles)
            for id in self.mutation["ID"].split(","):
                if id in self.seen_ids:
                    continue
                else:
                    self.seen_ids.add(id)
            if not self.mutation.split_info():  # If the information on the mutation type (intron, exon...) is missing
                continue
            if not self.mutation.rare_variant(threshold=MAF, populations_to_consider=self.valid_ancestries):
                continue
            if self.mutation.add_frequencies() and Genotypes.CompoundHeterozygotes in self.valid_genotypes:
                raise IOError("Cannot search for compound heterozygotes in unphased data")
            if self.check_exome is not None:
                if not self.check_exome.check_mutation(self.mutation):
                    continue
            self.get_coding_genes()
            for mutation in self.mutations:
                self.parse_new_mutation(mutation)
                self.remove_var_patients(mutation)
                if len(mutation["PATIENTS"]) > 0:
                    self.output_file.write(repr(mutation))

    def get_coding_genes(self):
        """
        Adds the "CODING_GENE" column.
        Each mutation can have multiple coding genes
        This function takes self.mutation, and finds all the coding genes for it
        It the creates copies of self.mutation, each one with one of the coding genes
        """
        self.mutations = []
        coding_genes = self.mutation.get_coding_genes(mutation_effects=self.valid_mutation_effects,
                                                      mutation_effect_lookup=self.mutation_effect_lookup,
                                                      search_level=self.search_level)
        for coding_gene in coding_genes:
            new_mutation = Mutation(self.mutation)
            new_mutation["CODING_GENE"] = coding_gene
            self.mutations.append(new_mutation)
        self.mutation = None

    def parse_new_mutation(self, mutation):
        """
        The goal of this function is to
        concatenates mutations that affect the same gene
        For example, if we have two mutations, that look like this
        ID      Coding_gene     Patient1    Patient2
        var1    gene1           1|0         1|1
        var2    gene1           0|1         0|0
        Our output would like this if Genotypes=[Genotypes.ComplexHeterozygous]
        Coding_gene     Patient1    Patient2
        gene1           1           0
        and this if Genotypes=[Genotypes.Homozygous]
        Coding_gene     Patient1    Patient2
        gene1           0           1

        This function modifies the self.ch_counts dictionary.
        If a mutation that affects this coding_gene has already been seen, it modifies that entry with the new
        information.  If not, it creates an entry.
        """
        coding_gene = mutation["CODING_GENE"]
        if coding_gene in self.ch_counts:
            self.ch_counts[coding_gene].parse_new_mutation(mutation)
        else:
            self.ch_counts[coding_gene] = CodingGeneMutation(self.ch_count_header, coding_gene, self.valid_genotypes)
            self.ch_counts[coding_gene].parse_new_mutation(mutation)

    def remove_var_patients(self, mutation):
        """
        This concatenates the patients affected by a mutation
        Beforehand, they were stored as
                patient1    patient2    patient3
        mut1    0|0         1|0         1|1
        Now, they are stored as
                PATIENTS
        mut1    patient2,patient3
        """
        patients_affected_by_mutation = []
        for patient in self.header.patient_columns:
            if mutation.data[patient][0] != "0" and mutation.data[patient][2] != "0":
                patients_affected_by_mutation.append(patient)
            del mutation.data[patient]
        mutation["PATIENTS"] = ",".join(patients_affected_by_mutation)

    def save_count_file(self, count_file_path):
        """
        Writes all count files (CodingGeneMutation, CodingGeneMutationHeader) to count_file_path
        """
        with open(count_file_path, 'w') as f:
            f.write(repr(self.ch_count_header))
            for coding_gene, ch_count in self.ch_counts.items():
                f.write(repr(ch_count))


def decode(s):
    encoding = "utf-8"
    if isinstance(s, bytes):
        return s.decode(encoding)
    else:
        return s


def read_parameter_file(parameter_file_path):
    valid_mutation_effects = []
    valid_ancestries = []
    valid_genotypes = []
    valid_search_levels = []
    MAF = None
    with open(parameter_file_path, 'r') as f:
        state = None
        enumeration = None
        valid_list = None
        for line in f:
            line = line.strip().upper()
            if line in ["MUTATION_EFFECTS", "ANCESTRY", "GENOTYPES", "SEARCH_LEVEL", "MINOR_ALLELE_FREQUENCY"]:
                state = line
            else:
                if state is None:
                    raise IOError("Needs header specifying Mutation_Effects, Ancestry, Genotypes, Search Level, or MAF")
                else:
                    if line.count(" ") == 1:
                        name = line.split(" ")[0]
                        if state == "MINOR_ALLELE_FREQUENCY":
                            if name == "MAF":
                                MAF = float(line.split(" ")[1])
                        if state == "MUTATION_EFFECTS":
                            enumeration = MutationEffects
                            valid_list = valid_mutation_effects
                        elif state == "ANCESTRY":
                            enumeration = Ancestry
                            valid_list = valid_ancestries
                        elif state == "GENOTYPES":
                            enumeration = Genotypes
                            valid_list = valid_genotypes
                        elif state == "SEARCH_LEVEL":
                            enumeration = SearchLevel
                            valid_list = valid_search_levels
                        for e_name, e_member in enumeration.__members__.items():
                            if e_name.lower() == name.lower():
                                valid_list.append(enumeration(e_member))
    if len(valid_search_levels) != 1:
        raise IOError("Search level in parameter file must be set to either Gene or Transcript")
    return valid_mutation_effects, valid_ancestries, valid_genotypes, valid_search_levels[0], MAF

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Search for phenotypic differences based on a vcf file")
    parser.add_argument('--input', '-i', help='Input vcf file', required=True)
    parser.add_argument('--output', '-o', help='Output vcf file', required=True)
    parser.add_argument('--count', '-c', help='Count output file', required=True)
    parser.add_argument('--parameters', '-p', help='Parameter file', required=True)
    parser.add_argument('--bed', '-b', help='Optional bed file')
    parser.add_argument('--subjects', '-s', help='Optional csv file with subject and population information')
    args = parser.parse_args()

    valid_mutation_effects, valid_ancestries, valid_genotypes, search_level, MAF = read_parameter_file(args.parameters)
    d = Data(annotated_file_path=args.input, output_file_path=args.output, count_file_path=args.count,
             subject_info_file_path=args.subjects, bed_file_path=args.bed,
             valid_mutation_effects=valid_mutation_effects, valid_ancestries=valid_ancestries,
             valid_genotypes=valid_genotypes, search_level=search_level, MAF=MAF)