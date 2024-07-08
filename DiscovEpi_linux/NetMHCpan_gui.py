#################################################
# Cedric Mahncke
# cedric.mahncke@leibniz-liv.de
# DiscovEpi
# Version 1.1
# Automatically retrieve protein data,
# predict corresponding epitopes and produce
# an epitope map for whole proteomes.
# 08.07.2024: Product.
#################################################

import re
import time
import operator
import requests as req
import utilities_gui as ut
from collections import Counter
import addtional_guis as add


# To fill data dict with given strings.
def fill_data(s, a, length):
    data = {
        'method': 'netmhcpan',
        'sequence_text': s,
        'allele': a,
        'length': length
    }
    return data


# Main function to execute prediction by NetMHCpan.
def exec_netmhcpan(uniprot, params, directory, progress, worker):
    start = time.time()
    progress.emit(1)
    # Dict containing sequences.
    seq_dict = uniprot[0]
    # Dict containing signal peptide ranges.
    sig_dict = uniprot[1]
    # Dict containing keys that have the same sequence.
    same_seq = uniprot[3]
    # Boolean about deletion of redundant sequences.
    redundant_del = uniprot[4]
    # Filename
    namestem = uniprot[5]
    # Dict with indices of epitopes located in the signal peptide.
    in_signal = {}
    # List containing all epitopes.
    epi_list = []
    # Storage for prediction data by location.
    epi_data = {}
    # Storage for whole prediction data.
    data = {}
    # Set parameters.
    allele, length, threshold = params
    params = {"Allel": allele, "Length": length, "Threshold": str(threshold)}

    # Calculate amount of predictions for providing progress time and percentage.
    # Prediction for redundant sequences does not affect progress time since their data is just copied.
    count_seqs = sum(len(seqs) for seqs in seq_dict.values())
    count_same_seqs = sum(len(seqs) for seqs in same_seq.values())

    if worker.isInterruptionRequested():
        return -3

    if not redundant_del:
        count_all_predictions = count_seqs - count_same_seqs
    else:
        count_all_predictions = count_seqs
    count = 0

    progress.emit(10)
    # Iterate over locations.
    for loc_key in seq_dict:

        if worker.isInterruptionRequested():
            return -3

        epi_data[loc_key] = {}
        in_signal[loc_key] = {}

        temp_keys = list(seq_dict[loc_key].keys())

        # Iterate over IDs of the proteins.
        i=0
        while temp_keys:

            if worker.isInterruptionRequested():
                return -3

            i+=1
            for gene_id in temp_keys:

                if worker.isInterruptionRequested():
                    return -3

                epi_data[loc_key][gene_id] = []
                in_signal[loc_key][gene_id] = []
                protein_name = seq_dict[loc_key][gene_id][1]

                # Check if an ID leading to an equal sequence exists.
                # If so, paste its prediction data to this entry.
                if gene_id in same_seq[loc_key]:
                    equal_seq_id = same_seq[loc_key][gene_id]
                    if equal_seq_id in data:
                        epi_list.extend(data[equal_seq_id][0])
                        epi_data[loc_key][gene_id] = data[equal_seq_id][1]
                        in_signal[loc_key][gene_id].extend(data[equal_seq_id][2])
                        count += 1
                        print(gene_id, "redundant sequence")
                        temp_keys.remove(gene_id)
                        continue

                # Specify the parameters and POST it for prediction.
                seq = seq_dict[loc_key][gene_id][0]
                pred_params = fill_data(seq, allele, length)
                header = {
                    'User-Agent': 'DiscovEpi Epitope Retrieval',
                    'From': 'cedric.mahncke@leibniz-liv.de'
                }

                try:
                    response = req.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/', data=pred_params,
                                        headers=header)
                except:
                    print("Error at: " + gene_id)
                    continue

                if response.ok:
                    # Filter important data, calculate densities and add everything to storages.
                    ut.request_handle(method="NetMHCpan", request=response, indices=[loc_key, gene_id],
                                      pred_vars=[len(seq), length, threshold, allele],
                                      epi_storages=[epi_list, epi_data, data], signal_storages=[sig_dict, in_signal],
                                      protein_name=protein_name)
                    temp_keys.remove(gene_id)
                    count += 1
                else:
                    # Print note at console and output file if the request is not okay.
                    print(loc_key, gene_id, "\tBad Answer, waiting 5 seconds and trying again.")
                    time.sleep(5)
                    continue

                progress.emit(10+88*count/count_all_predictions)

    # Check if epi_list is empty.
    if not epi_list:
        print("No epitopes found.")
        return -1

    # Count appearance of each epitope.
    frequencies = dict(sorted(Counter(epi_list).items(), key=operator.itemgetter(1), reverse=True))

    # Filename
    allele_sub = re.sub(r"[-\*:]", "_", allele)
    filename = directory + "\\" + re.sub("_+", "_", "_".join(["nmp", namestem, loc_key, allele_sub])) + ".xlsx"
    print(filename)
    header = 'id,protein,score,density,sig density,position; epitope; score'
    ut.print_sheet(header, params, epi_data, frequencies, filename, "NetMHCPan-4.0", in_signal)

    end = time.time()
    print("Time elapsed in NetMHCpan: ", add.calc_time(end - start))
    return epi_data, filename
