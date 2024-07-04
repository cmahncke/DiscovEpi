#################################################
# Cedric Mahncke
# s-cemahn@uni-greifswald.de
# DiscovEpi
# Version 1.0
# A tool to automatically retrieve protein data,
# predict corresponding epitopes and produce
# potential epitope binding maps for whole proteome.
# 20.11.2023: Product.
#################################################


import requests as req
import addtional_guis as add
import utilities_gui as ut
from collections import Counter
import operator
import re
import time


# To make the request query by given arguments.
def make_query_args(params):
    try:
        # For every location a query is made and stored with the same location as key.
        queries = {}

        #remove slash from organism name
        params[0] = params[0].replace("/", "")
        if "," in params[1]:
            locs = params[1].split(",")
        elif params[1] == "":
            locs = ["All"]
        else:
            locs = [params[1]]

        if params[2]:
            params[2] = " AND reviewed:true"
        else:
            params[2] = ""

        for loc_key in locs:
            if loc_key == "All":
                queries[loc_key] = 'organism_name:' + params[0] + params[2]
            else:
                queries[loc_key] = 'organism_name:' + params[0] + ' AND cc_scl_term:' + loc_key + params[2]

        params = {'Organism': params[0], 'Locations': params[1], 'Reviewed': params[2].replace("AND reviewed:", "")}
        return queries, params
    except ValueError:
        print("Value error occurred when making query! ...try again.")


def progress_changed(progress, value):
    progress.setValue(value)


# Based on the entries made, the function executes the request from uniprot.org, sorts out entries for different
# organism parameter and stores whole data indexed by location and ID and also the ID of an entry with the same sequence
# if there is one indexed by the current ID.
def search_uniprot(queries, columns, params):
    try:
        header = {
            'User-Agent': 'DiscovEpi Protein Retrieval',
            'From': 'cedric.mahncke@leibniz-liv.de'
        }
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
                       'format': 'tsv',
                       'fields': columns}
            try:
                print(payload)
                request = req.get(base + resource, params=payload, headers=header)
            except:
                raise ConnectionError

            print(request.url)
            result[loc_key] = {}
            same_seq[loc_key] = {}
            if request.ok:
                # In tab format the data is stored in a nested dictionary with the first column, id in
                # example, as secondary key.
                # Delete header and tail of the request.
                entries = request.text.splitlines()[1:]
                print(len(entries))
                if len(entries) == 0:
                    raise TypeError
                for entry in entries:
                    split = entry.split("\t")
                    # Filter not wanted strains
                    strain = split[-3].split(",")[-1]
                    strain = strain.replace(' / ', '/')
                    if params['Organism'] not in strain:
                        print(f'{split[0]} skipped because strain {strain} is not {params["Organism"]}!')
                        continue
                    if params['Organism'] == 'SARS-CoV' and 'SARS-CoV-2' in strain:
                        print(f'{split[0]} skipped because strain {strain} is not {params["Organism"]}!')
                        continue
                    # Store keys with equal sequences
                    seq = split[-1]
                    gene_id = split[0]
                    if seq in sequences:
                        same_seq[loc_key][gene_id] = seq_key[seq]
                    else:
                        sequences.append(seq)
                        seq_key[seq] = gene_id

                    # Strip unnecessary text from location
                    split[-4] = split[-4].strip("SUBCELLULAR LOCATION: ")

                    # Filter signal range
                    signal_evidence = split[-2].split(";")
                    if signal_evidence != ['']:
                        signal = signal_evidence[0].split()[1]
                        signal_evidence[0] = signal
                        split[-2] = signal.split(".")[-1]
                    else:
                        split[-2] = "N/A"

                    # Nested dictionary with locations as primary keys and the first column, e.g. id, as secondary.
                    result[loc_key][gene_id] = split[1:]
            else:
                raise TypeError
        if result.values() == {}:
            raise TypeError

        return result, same_seq
    except ConnectionError:
        print("Connection error occurred in UniProt request function! ...try again.")
        raise ConnectionError


# Execute the data collection and preparing for proteome output
def exec_uniprot(parameter, del_redundant, directory, progress, worker):
    start = time.time()
    # Universal algorithm with parameters input by user.
    progress.emit(1)

    q, params = make_query_args(parameter)
    progress.emit(20)

    columns = "id,protein_name,gene_names,length,cc_subcellular_location,organism_name,ft_signal,sequence"
    columns_txt = "id,protein,gene,length,location,organism,signal_sequence,sequence"

    if worker.isInterruptionRequested():
        return -3

    try:
        data, same_seq = search_uniprot(q, columns, params)
    except TypeError:
        progress.emit(0)
        print("UniProt Search Error")
        return -1
    except ConnectionError:
        progress.emit(0)
        print("UniProt Connection Error")
        return -2

    # Storages for passing sequences and signal peptides to prediction functions, counting the sequences and
    # provide all deleted sequences.
    seq_dict = {}
    seq_list = []
    sig_dict = {}
    popped_seqs = {}
    progress.emit(40)

    # Based on the requested uniprot data the storages are filled. First keyed by location and then by ID.
    for loc_key in data:

        if worker.isInterruptionRequested():
            return -3

        if del_redundant:
            popped_seqs[loc_key] = ut.del_redundant(data[loc_key], same_seq[loc_key])

        seq_dict[loc_key] = {}
        sig_dict[loc_key] = {}

        for gene_id in data[loc_key]:

            if worker.isInterruptionRequested():
                return -3


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
    ut.print_sheet(columns_txt, params, data, frequencies, filename, "UniProt")
    progress.emit(99)

    end = time.time()
    print("Time elapsed in UniProt: ", add.calc_time(end - start))
    return seq_dict, sig_dict, popped_seqs, same_seq, del_redundant, namestem
