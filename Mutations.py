from Enumerations import Ancestry, MutationEffects, Genotypes, SearchLevel, SuperPopulations
import os
import re
from collections import defaultdict


class Header:
    def __init__(self, *args):
        if len(args) == 1:
            if isinstance(args[0], str):
                """
                This option reads the column headers in from an input file
                input_titles: read in from vcf file
                patient_columns: columns that match patient ids
                output_titles: to be output to a vcf file
                """
                header_line = args[0]
                self.input_titles = header_line.split("\t")
                info_columns = ["AA", "AC", "AF", "AFR_AF", "AMR_AF", "AN", "DP", "EAS_AF", "EUR_AF", "HRun", "NS",
                                "SAS_AF", "VT", "EFF"]
                self.patient_columns = [x for x in self.input_titles if bool(re.search(r'\d', x))]
                self.output_titles = self.input_titles[0:1] + ["PATIENTS"] + self.input_titles[1:4] + ["CODING_GENE"]\
                    + self.input_titles[4:7] + info_columns + \
                    [x for x in self.input_titles[7:] if x not in self.patient_columns] + \
                    ["AC", "AMax", "AF"]
                #self.var_columns = [x for x in self.output_titles if x not in self.count_columns]
            elif isinstance(args[0], Header):
                """
                This option creates a copy of an instance of the Header class
                """
                header = args[0]
                self.input_titles = header.input_titles
                self.output_titles = header.output_titles
                self.patient_columns = header.patient_columns
                #self.var_columns = header.var_columns
            else:
                raise ValueError("Unrecognized arguments")
        else:
            raise ValueError("Wrong number of arguments")

    def __getitem__(self, item):
        return self.output_titles[item]

    def __len__(self):
        return len(self.output_titles)

    def __repr__(self):
        return "\t".join(self.output_titles) + os.linesep


class Mutation:
    def __init__(self, *args):
        if len(args) == 3:
            """
            This option loads a line of the vcf file as a mutation
            The header (titles) have been put into self.header
            One line, detailing one mutation, will be put into this object, and placed in Data.mutation
            """
            header, line, abbreviated_titles = args
            if isinstance(header, Header) and isinstance(line, str) and isinstance(abbreviated_titles, dict):
                columns = line.split("\t")
                self.data = {}
                for title, column in zip(header.input_titles, columns):
                    self[title] = column
                self.header = header
                self.abbreviated_titles = abbreviated_titles
                try:
                    id = self["ID"]
                except KeyError as e:
                    print(line)
                    raise e
            else:
                raise ValueError("Unrecognized arguments")
        elif len(args) == 1:
            """
            This option create a copy of an instance of the Mutation class
            I couldn't figure out how to make __copy__ work
            """
            mutation = args[0]
            if isinstance(mutation, Mutation):
                self.header = mutation.header
                self.abbreviated_titles = mutation.abbreviated_titles
                self.data = dict(mutation.data)
            else:
                raise ValueError("Unrecognized arguments")
        else:
            raise ValueError("Wrong number of arguments")

    def add_frequencies(self):
        """
        This function computes some simple statistics on each mutation.
        Computes the number of mutated alleles, the number of wild-type alleles, and the frequency of mutations

        Sometimes the reference and alternate allele are flipped.  To counteract this, if more than 50% of the
        alleles are alternate alleles, this function will modify the data so that the alternate allele
        is the reference allele

        Some patients have more than one alternate allele.
        If so, the allele count will be the frequency of all alternate alleles
        """
        unphased_found = False
        allele_counts = defaultdict(int)
        for patient in self.header.patient_columns:
            if self[patient][1] == "/":  # If it's unphased, just ignore it, and assume both alleles are reference
                unphased_found = True
            allele_counts[self[patient][0]] += 1
            allele_counts[self[patient][-1]] += 1
        self['AMax'] = 2 * len(self.header.patient_columns)
        self['AC'] = self['AMax'] - allele_counts[0]

        reference_allele_number = -1  # This code chunk identifies which allele has the greatest frequency
        reference_count = -1
        for allele_number, allele_count in allele_counts.items():
            if allele_count > reference_count and allele_number != ".":
                reference_allele_number = int(allele_number)
                reference_count = allele_count

        if reference_allele_number != 0:  # If the reference allele is not the most common
            alternate_alleles = self["ALT"].split(",")
            self["REF"], alternate_alleles[reference_allele_number - 1] = \
                alternate_alleles[reference_allele_number - 1], self["REF"]
            self["ALT"] = ",".join(alternate_alleles)
            self['AC'] = self['AMax'] - allele_counts[reference_allele_number]
            reference_allele_number = str(reference_allele_number)
            for patient in self.header.patient_columns:
                new_genotype = ""
                allele = self[patient][0]
                if allele == '0':
                    new_genotype += reference_allele_number
                elif allele == reference_allele_number:
                    new_genotype += "0"
                else:
                    new_genotype += allele
                new_genotype += "|"
                allele = self[patient][2]
                if allele == '0':
                    new_genotype += reference_allele_number
                elif allele == reference_allele_number:
                    new_genotype += "0"
                else:
                    new_genotype += allele
                self[patient] = new_genotype
        self['AF'] = float(self['AC']) / float(self['AMax'])
        return unphased_found

    def split_info(self):
        """
        Column "INFO" looks like "AF=0.2;EUR_AF=0.3;EAS_AF=0.1..."
        split this into
        Column "AF" = 0.2
        Column "EUR_AF" = 0.3
        Column "EAS_AF" = 0.1
        ...
        Returns whether it created the columns "ANN" or "EFF",
        which are used to identify the effects of the mutation on the gene
        """
        info = self["INFO"].split(";")
        del self.data["INFO"]
        for datum in info:
            datum = datum.strip('|').split("=")
            if len(datum) == 2:
                self[datum[0]] = datum[1]
        if "EFF" not in self and "ANN" not in self:
            print("Couldn't determine effect of mutation for %s" % self["ID"])
            return False
        else:
            return True

    def rare_variant(self, threshold, populations_to_consider=Ancestry.all()):
        """
        This function returns whether or not this mutation is rare;
        if it has an incidence of MAF or greater in the populations
        However, mutation.rare_variant takes an optional argument which can specify which populations should be checked
        :populations_to_consider: A list of objects of type Ancestry
        :return: Boolean, if this mutation is "rare"; <.02 for all populations
        """
        try:
            if Ancestry.Overall in populations_to_consider:
                if not self.rare_variant_population(threshold, 'AF'):
                    return False
            if Ancestry.African in populations_to_consider:
                if not self.rare_variant_population(threshold, 'AFR_AF'):
                    return False
            if Ancestry.American in populations_to_consider:
                if not self.rare_variant_population(threshold, 'AMR_AF'):
                    return False
            if Ancestry.EastAsian in populations_to_consider:
                if not self.rare_variant_population(threshold, 'EAS_AF'):
                    return False
            if Ancestry.European in populations_to_consider:
                if not self.rare_variant_population(threshold, 'EUR_AF'):
                    return False
            if Ancestry.SouthAsian in populations_to_consider:
                if not self.rare_variant_population(threshold, 'SAS_AF'):
                    return False
            if Ancestry.ThousandGenomesOverall in populations_to_consider:
                if not self.rare_variant_population(threshold, 'dbNSFP_1000Gp1_AF'):
                    return False
            if Ancestry.ThousandGenomesAfrican in populations_to_consider:
                if not self.rare_variant_population(threshold, 'dbNSFP_1000Gp1_AFR_AF'):
                    return False
            if Ancestry.ThousandGenomesAmerican in populations_to_consider:
                if not self.rare_variant_population(threshold, 'dbNSFP_1000Gp1_AMR_AF'):
                    return False
            if Ancestry.ThousandGenomesAsian in populations_to_consider:
                if not self.rare_variant_population(threshold, 'dbNSFP_1000Gp1_ASN_AF'):
                    return False
            if Ancestry.ThousandGenomesEuropean in populations_to_consider:
                if not self.rare_variant_population(threshold, 'dbNSFP_1000Gp1_EUR_AF'):
                    return False
        except ValueError:  # If one of the allele frequencies cannot be cast to a float
            return False
        return True

    def rare_variant_population(self, threshold, population):
        if population not in self:
            return True
        return not (threshold < float(self[population]) < (1 - threshold))

    def get_coding_genes(self, mutation_effects=None, mutation_effect_lookup=None, search_level=SearchLevel.Gene):
        """
        This code is very annotation dependent, and may have to be changed depending on how you annotated your files

        If search_level is set to SearchLevel.Transcript, then this will deal with the transcript, not the gene

        Each mutation has a column, labeled "EFF" or "ANN" detailing how this mutation will effect the genes
        This column contains a lot of information, but we only care about the affected gene and how it is affected
        If the column is called "EFF", it was created by an old version of the annotation software, and it'll look like:
            ______              __________
            INTRON(MODIFIER|||||AP000525.8|processed_transcript|NON_CODING|ENST00000413768|7)
            gene effect         gene                                       transcript
        If the column is called "ANN", it was annotated by a newer software, and it'll look like:
              __________________          _______________
            A|SYNONYMOUS_VARIANT|LOW|TPTE|ENSG00000166157|TRANSCRIPT|ENST00000342420|PROTEIN_CODING|16/22|C.939C>T|P.ALA313ALA|1269/2036|939/1542|313/513||
              gene effect                 gene                       transcript
        This function only returns coding_genes that will by changed by the mutation
        As such, the user can specify if they're looking for mutations that are changing Exons, Introns, Synonymous...
            (Further options in Enumerations.MutationEffects)
        This option, to determine which mutation effects are acceptable, is passed in as mutation_effects

        This function also allows users to specify a mutation_effect_lookup, to get around the problem that sometimes
        the same mutation effect is called different names.
        For example, sometimes Downstream is called "DOWNSTREAM" and sometimes "DOWNSTREAM_GENE_VARIANT"
        Pass in a dictionary with two entries;
            d["DOWNSTREAM"] = MutationEffects.DOWNSTREAM
            d["DOWNSTREAM_GENE_VARIANT"] = MutationEffects.DOWNSTREAM
        If no dictionary is passed in, it will default to MutationEffects.str_lookup()
        """
        if mutation_effects is None:
            mutation_effects = MutationEffects.default()
        if mutation_effect_lookup is None:
            mutation_effect_lookup = MutationEffects.str_lookup()
        coding_genes = set()
        if "EFF" in self:
            effects = self["EFF"].split(",")
            for effect in effects:
                effect = effect.upper()
                gene_effect = effect.split("(")[0]
                if search_level == SearchLevel.Gene:
                    gene = effect.split("|")[5]
                else:
                    gene = effect.split("|")[8]
                for each_effect in gene_effect.split("&"):
                    if each_effect not in mutation_effect_lookup:
                        each_effect = each_effect.split("+")[0]
                    if mutation_effect_lookup[each_effect] in mutation_effects:
                        coding_genes.add(gene)
        elif "ANN" in self:
            effects = self["ANN"].split(",")
            for effect in effects:
                effect = effect.upper().split("|")
                gene_effect = effect[1]
                if search_level == SearchLevel.Gene:
                    gene = effect[3]
                else:
                    gene = effect[6]
                for each_effect in gene_effect.split("&"):
                    if each_effect not in mutation_effect_lookup:
                        each_effect = each_effect.split("+")[0]
                    if each_effect in mutation_effect_lookup:
                        if mutation_effect_lookup[each_effect] in mutation_effects:
                            coding_genes.add(gene)
                    else:
                        print(effect)
        return coding_genes

    def __repr__(self):
        return '\t'.join([str(self.data[title]) if title in self.data.keys() else ""
                          for title in self.header.output_titles]) + os.linesep

    def __getitem__(self, item):
        """
        This takes advantage of the complicated abbreviated title lookup.
        If you're searching for item...
        -Check if item is in the dictionary
            -If so, return it
            -If not, check if item is in the abbreviated lookup dictionary (ie, 'Ancestry' for 'Ancestry Alleles')
                -If so, lookup the full key, then use that to lookup the value
                -If not, look to see if item is a new, previously unseen, abbreviation
                    -If so, put a new entry in the abbreviated lookup dictionary
                        (ie self.abbreviated_titles['Ancestry'] = 'Ancestry Alleles')
                    -If item matches too many or too few entries, raise an Exception
        """
        if item in self.data:
            return self.data[item]
        else:
            if item in self.abbreviated_titles:
                return self.data[self.abbreviated_titles[item]]
            else:
                possible_key_matches = []
                for key in self.data.keys():
                    if item in key:
                        possible_key_matches.append(key)
                if len(possible_key_matches) == 1:  # This is a substring of one and only one column
                    self.abbreviated_titles[item] = possible_key_matches[0]
                    return self.data[possible_key_matches[0]]
                elif len(possible_key_matches) == 0:
                    raise KeyError("No matches found.  Misspelling of %s?" % item)
                else:
                    raise KeyError("Matches more than one entry")

    def __setitem__(self, key, value):
        self.data[key] = value

    def __contains__(self, item):
        return item in self.data


class PatientGenotype:
    """
    This class takes care of one coding_gene for one patient, for one list of genotypes
    For example, AP00025.1 for HG01302 and [ComplexHeterozygotes]
    Every time you see a new mutation that will affect AP00025.1, call parse_new_alleles
    self.value holds whether or not this coding_gene for this patient is one of the types in valid_genotypes
    """
    def __init__(self, valid_genotypes):
        self.value = False  # Represents whether or not this patient has this genotype
        self.maternal_seen = False
        self.paternal_seen = False
        self.valid_genotypes = valid_genotypes

    def parse_new_alleles(self, alleles):
        """
        Update self.value with new information
        """
        if self.value is False:
            maternal_allele = alleles[0] != "0"
            paternal_allele = alleles[2] != "0"  # False is wildtype, True is mutated
            if Genotypes.OneMutation in self.valid_genotypes:
                if maternal_allele or paternal_allele:
                    self.value = Genotypes.OneMutation
                    return
            if Genotypes.TwoMutations in self.valid_genotypes:
                if maternal_allele and (paternal_allele or self.paternal_seen):
                    self.value = Genotypes.TwoMutations
                    return
                if paternal_allele and (maternal_allele or self.maternal_seen):
                    self.value = Genotypes.TwoMutations
                    return
            if Genotypes.Homozygous in self.valid_genotypes:
                if maternal_allele and paternal_allele:
                    self.value = Genotypes.Homozygous
                    return
            if Genotypes.CompoundHeterozygotes in self.valid_genotypes:
                if maternal_allele and self.paternal_seen:
                    self.value = Genotypes.CompoundHeterozygotes
                    return
                if paternal_allele and self.maternal_seen:
                    self.value = Genotypes.CompoundHeterozygotes
                    return
            self.maternal_seen = self.maternal_seen or maternal_allele
            self.paternal_seen = self.paternal_seen or paternal_allele

    def __repr__(self):
        if self.value is False:
            return "0"
        elif self.value is Genotypes.OneMutation:
            return "o"
        elif self.value is Genotypes.TwoMutations:
            return "t"
        elif self.value is Genotypes.CompoundHeterozygotes:
            return "c"
        elif self.value is Genotypes.Homozygous:
            return "h"


class CodingGeneMutation(Mutation):
    """
    This class holds genotype changes for every patient for one coding gene
    This class will have one PatientGenotype object for every patient,
    These PatientGenotype objects will be updated with every new mutation affecting the given coding_gene
    """
    def __init__(self, ch_count_header, coding_gene, valid_genotypes):
        self.header = ch_count_header
        self.abbreviated_titles = {}
        self.valid_genotypes = valid_genotypes

        self.data = {}
        for key in self.header.count_columns:
            self.data[key] = 0
        self.data['CODING_GENE'] = coding_gene
        for patient in self.header.patient_columns:
            self.data[patient] = PatientGenotype(self.valid_genotypes)

    def parse_new_mutation(self, mutation):
        for patient in self.header.patient_columns:
            self.data[patient].parse_new_alleles(mutation[patient])

    def __repr__(self):
        return '\t'.join(
            [str(self.data[title]) if title in self.data.keys() else "" for title in self.header.output_titles])\
               + os.linesep


class CodingGeneMutationHeader(Header):
    """
    Creates a header for the Counts
    """
    def __init__(self, header, my_patients):
        self.count_columns = ["CODING_GENE", "COUNT"]
        for key, values in my_patients.possibilities.items():
            for value in values:
                self.count_columns += [key + "=" + value]

        self.output_titles = self.count_columns + header.patient_columns
        self.patient_columns = header.patient_columns
        self.my_patients = my_patients

    def __repr__(self):
        patient_header = '\t'.join(self.output_titles)
        subject_info_headers = []
        for key in self.my_patients.possibilities:
            s = '\t' * len(self.count_columns)
            s += '\t'.join([self.my_patients.patients[patient].data[key] for patient in self.patient_columns])
            subject_info_headers.append(s)
        return patient_header + os.linesep.join(subject_info_headers)
