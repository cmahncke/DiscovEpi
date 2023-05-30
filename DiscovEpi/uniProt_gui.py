#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscoEpi
# A Tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 05.05.2023: Product.
#################################################


import requests as req
import utilities_gui as ut
from collections import Counter
import operator
import re


# To make the request query by given arguments.
def make_query_args(params):
    try:
        # For every location a query is made and stored with the same location as key.
        queries = {}
        if "," in params[1]: locs = params[1].split(",")
        elif params[1] == "": locs = ["All"]
        else: locs = [params[1]]

        if params[2]: params[2] = "true"
        else: params[2] = "false"

        for loc_key in locs:
            if loc_key == "All":
                queries[loc_key] = 'organism_name:' + params[0] + ' AND reviewed:' + params[2]
            else:
                queries[loc_key] = 'organism_name:' + params[0] + ' AND cc_scl_term:' + loc_key + ' AND reviewed:' + params[2]

        params = {'Organism': params[0], 'Locations': params[1], 'Reviewed': params[2]}
        return queries, params
    except ValueError:
        print("Value error occurred when making query! ...try again.")


def progress_changed(progress, value):
    progress.setValue(value)


# Based on the entries made, the function executes the request from uniprot.org, sorts out entries for different
# organism parameter and stores whole data indexed by location and ID and also the ID of an entry with the same sequence
# if there is one indexed by the current ID.
def search_uniprot(queries, columns, form, params):
    try:
        base = 'https://rest.uniprot.org'
        resource = '/uniprotkb/stream?'
        result = {}
        same_seq = {}

        # For every location there is a request loaded into a dictionary keyed by locations.
        for loc_key in queries:

            # For storing keys to redundant sequences
            sequences = []
            seq_key = {}
            payload = {'query': queries[loc_key],
                       'format': form,
                       'fields': columns}
            request = req.get(base + resource, params=payload)
            print(request.text)
            print(request.url)
            result[loc_key] = {}
            same_seq[loc_key] = {}
            if request.ok:
                if form == "tsv":
                    # In tab format the data is stored in a nested dictionary with the first column, id in
                    # example, as secondary key.
                    # Delete header and tail of the request.
                    entries = request.text.splitlines()[1:]
                    print(entries)
                    for entry in entries:
                        split = entry.split("\t")
                        # Filter not wanted strains
                        strain = split[-3].split(",")[-1]
                        con = False
                        for word in params['Organism'].split():
                            if word not in strain:
                                print(word, "not in strain: ", strain, "; continue", split[0])
                                con = True
                        if con: continue
                        # Store keys with equal sequences
                        seq = split[-1]
                        gene_id = split[0]
                        if seq in sequences:
                            same_seq[loc_key][gene_id] = seq_key[seq]
                        else:
                            sequences.append(seq)
                            seq_key[seq] = gene_id

                        # Filter signal range
                        signal_evidence = split[-2].split(";")
                        if signal_evidence != ['']:
                            signal = signal_evidence[0].split()[1]
                            signal_evidence[0] = signal
                            split[-2] = signal.split(".")[-1]
                        else: split[-2] = "N/A"

                        # Nested dictionary with locations as primary keys and the first column, e.g. id, as secondary.
                        result[loc_key][gene_id] = split[1:]

                else:
                    result[loc_key] = request.text

            else:
                print(request.text)
                print("Error: Request not okay!")
                return -1

        return result, same_seq
    except ConnectionError:
        print("Connection error occurred in UniProt request function! ...try again.")


# Execute the data collection and preparing for proteome output
def exec_uniprot(parameter, out_format, del_redundant, directory, progress, QMainWindow):
    # Universal algorithm with parameters input by user.
    progress.emit(0)

    q, params = make_query_args(parameter)
    progress.emit(20)

    form = out_format.lower()
    if form == "tsv":
        columns = "id,protein_name,gene_names,length,cc_subcellular_location,organism_name,ft_signal,sequence"
        try:
            data, same_seq = search_uniprot(q, columns, form, params)
        except TypeError:
            progress.emit(0)
            print("Search Error")
            return -1

        # Storages for passing sequences and signal peptides to prediction functions, counting the sequences and
        # provide all deleted sequences.
        seq_dict = {}
        seq_list = []
        sig_dict = {}
        popped_seqs = {}
        progress.emit(40)

        # Based on the requested uniprot data the storages are filled. First keyed by location and then by ID.
        for loc_key in data:

            if del_redundant:
                popped_seqs[loc_key] = ut.del_redundant(data[loc_key], same_seq[loc_key])

            seq_dict[loc_key] = {}
            sig_dict[loc_key] = {}

            for gene_id in data[loc_key]:
                # Fill sequence dict with sequence and protein name.
                seq_dict[loc_key][gene_id] = [data[loc_key][gene_id][-1], data[loc_key][gene_id][0]]

                # Fill signal dict with signal range.
                sig_dict[loc_key][gene_id] = data[loc_key][gene_id][-2]

                # Fill sequence list with sequence.
                seq_list.append(data[loc_key][gene_id][-1])
        progress.emit(60)

        # The dict contains the number of appearances sorted decreasingly keyed by sequence.
        info = {"Warning:": "If redundant have been deleted only the first appearence of each sequence is stored."}
        frequencies = {**info, **dict(sorted(Counter(seq_list).items(), key=operator.itemgetter(1), reverse=True))}
        progress.emit(80)
        # Filename
        namestem = re.sub(r"[\s,.:\-/\\()]", "_", parameter[0])
        #if parameter[1]:
        #    namestem = "_".join([namestem,
        #                         re.sub(r"[\s,.:\-/\\()]", "_", parameter[1])])
        filename = directory + "\\" + re.sub("_+", "_", "_".join(["unp", namestem, loc_key])) + ".xlsx"
        # Function to print data in a spreadsheet is defined in 'utilities.py'.
        ut.print_sheet(columns, params, data, frequencies, filename, "UniProt")
        progress.emit(99)
        print(seq_dict)
        return seq_dict, sig_dict, popped_seqs, same_seq, del_redundant, namestem

    # If requested format is fasta columns are not used and predictions not made.
    elif form == "fasta":
        fasta_data = search_uniprot(queries=q, columns="", form=form, params=params)

        # For every location there is a .fsa file containing FASTA data to every protein.
        for loc_key in fasta_data:
            filename = "uniProt_" + loc_key + ".fsa"
            ut.print_file(fasta_data[loc_key], filename)

        return 0

    else:
        print("Error: Unknown file format!\nOnly 'tsv' and 'fasta' until now.")
        return -1
