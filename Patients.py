import re


class Patient:
    """
    Used to represent information about each patient
    Gender, ancestry....
    """
    def __init__(self, patient_name, titles=None, line=None):
        """
        Insert data into dictionary
        patient_name = 'HG0096'
        titles = 'Patient_name', 'Gender', 'Population'
        line = 'HG0096, M, GBR'

        causes
        self.patient_name = 'HG0096'
        self.data = ['Gender': 'M', 'Population': 'GBR']
        """
        self.data = {}
        self.patient_name = patient_name
        if titles is not None and line is not None:
            fields = line.strip().split(",")
            if len(titles) != len(fields):
                raise TypeError("Different title and line length.  "
                                "Please remove any extra commas from Subject Info file")
            for title, field in zip(titles, fields):
                self.data[title] = field
            for key in list(self.data.keys()):
                if self.data[key] == self.patient_name:
                    del self.data[key]

    def __hash__(self):
        return hash(self.patient_name)

    def __eq__(self, other):
        return self.patient_name == other.patient_name


class Patients:
    def __init__(self, patients=None, subject_info_file_path=None):
        """
        Creates a list of Patients from an input file or from an existing list
        Call this only with patients OR subject_info_file_path
        """
        self.possibilities = {} # Holds the possible entries for every header
        if patients is not None:
            self.patients = {}
            for patient_name in patients:
                self.patients[patient_name] = Patient(patient_name=patient_name)
        if subject_info_file_path is not None:
            self.patients = {}
            with open(subject_info_file_path) as subject_info:
                titles = subject_info.readline().strip().split(",")  # Get rid of header line
                for line in subject_info:
                    if bool(re.search(r'\d', line)):
                        p = Patient(patient_name=line.split(",")[0], titles=titles, line=line)
                        for key, value in p.data.items():
                            if key in self.possibilities:
                                if value not in self.possibilities[key]:
                                    self.possibilities[key].append(value)
                            else:
                                self.possibilities[key] = [value]
                        self.patients[p.patient_name] = p

    def __getitem__(self, item):
        return self.patients[item]

    def __contains__(self, item):
        return item in self.patients
