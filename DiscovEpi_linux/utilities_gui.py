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

import xlsxwriter
import numpy as np
import re
from datetime import datetime


# This is to filter important data from the requests based on the algorithm used. Also densities are being calculated
# and data added to storages inside the functions.
def request_handle(method, request, indices, pred_vars, epi_storages, signal_storages, protein_name):
    loc_key, gene_id = indices
    seq_len, length, threshold, allele = pred_vars
    epi_list, epi_data, global_data = epi_storages
    signal_dict, inside_signal = signal_storages

    temp_epi_list = []
    temp_score = {}
    temp_epi_data = {}

    lines = request.text.splitlines()[1:]
    for line in lines:
        # Filter important information.
        pred_data = line.split("\t")
        if len(pred_data) >= 10:
            position, epitope, score = pred_data[2], pred_data[5], pred_data[9]
        else:
            print("Data not complete for epitop with index", lines.index(line), "of protein", gene_id)
            break

        # Only take top scoring epitopes.
        if float(score) <= float(threshold):
            # Put the score in relation to the threshold and make it a relative scale.
            score = str(round((float(threshold) - float(score)) / float(threshold), 2))

            # Store index of epitopes which are located inside the signal peptide.
            if signal_dict[loc_key][gene_id].isnumeric() and int(position) <= int(signal_dict[loc_key][gene_id]):
                inside_signal[loc_key][gene_id].append(lines.index(line) + 3)

            # Add filtered information to list and dictionary about epitopes.
            temp_epi_list.append(epitope)
            temp_score[epitope] = float(score)
            temp_epi_data[epitope] = "; ".join([position, epitope, score])

    # Append protein name for information.
    epi_data[loc_key][gene_id].append(protein_name)

    # Calculate epitope density.
    # Overall.
    epitopes_bound = len(temp_epi_data)
    epitopes_total = seq_len - int(length) + 1

    overall_score = sum(temp_score.values())
    if epitopes_bound > 0:
        #relative_score = overall_score / epitopes_bound
        relative_score = overall_score / epitopes_total
        density = epitopes_bound / epitopes_total
    else:
        density = 0
        relative_score = 0
    # Store density with epitope data.
    epi_data[loc_key][gene_id].extend([relative_score, density])

    # Inside signal peptide.
    # If signal is given as number, not 'N/A'.
    if signal_dict[loc_key][gene_id].isnumeric():
        signal_length = int(signal_dict[loc_key][gene_id])
        signal_epitopes_bound = len(inside_signal[loc_key][gene_id])
        signal_density = signal_epitopes_bound / signal_length

        # Store signal density with epitope data.
        epi_data[loc_key][gene_id].append(signal_density)
    else:
        epi_data[loc_key][gene_id].append("N/A")

    epi_data[loc_key][gene_id].extend(temp_epi_data.values())
    epi_list.extend(temp_epi_list)
    global_data[gene_id] = [temp_epi_list, epi_data[loc_key][gene_id], inside_signal[loc_key][gene_id]]
    return 0


# Print text to file.
def print_file(text, file):
    with open(file, 'w') as fi:
        fi.write(text)
    fi.close()
    return 0


# Extend content of existing file with text.
def extend_file(text, file):
    with open(file, 'a') as f:
        f.write(text)
    f.close()
    return 0


# To write requested data into .xlsx files.
def print_sheet(header, params, str_dict, freq_dict, filename, method, *signal_mark):
    if signal_mark:
        in_signal = signal_mark[0]
    else:
        in_signal = None

    # Create workbook and add first sheet for overview.
    wb = xlsxwriter.Workbook(filename)
    sheet = wb.add_worksheet("method")
    form = wb.add_format()

    # Write first page about the method used, parameters set and frequencies of sequences / epitopes.
    form.set_bold()
    # Method.
    sheet.write(0, 0, method, form)
    # Parameters that has been set.
    sheet.write(2, 0, "Parameter:", form)
    for param in params:
        col_index = list(params.keys()).index(param)
        sheet.write(1, col_index + 1, param, form)
        if params[param].isnumeric():
            sheet.write(2, col_index + 1, int(params[param]), form.set_num_format("0"))
        else:
            sheet.write(2, col_index + 1, params[param])
    # Logs.
    sheet.write(3, 0, "Logs:", form)
    sheet.write(3, 1, datetime.now().date().strftime("%Y-%m-%d"))
    sheet.write(3, 2, datetime.now().time().strftime("%H:%M:%S"))
    # Frequencies.
    sheet.write(4, 0, "Frequencies over all loci:", form)
    for gene_id in freq_dict:
        row_index = list(freq_dict.keys()).index(gene_id)
        sheet.write(row_index + 5, 0, gene_id)
        if isinstance(freq_dict[gene_id], int):
            sheet.write(row_index + 5, 1, int(freq_dict[gene_id]))
        else:
            sheet.write(row_index + 5, 1, freq_dict[gene_id])

    # Link base.
    base = "http://www.uniprot.org/uniprot/"
    header = header.split(",")
    header.insert(-1, "link")

    # For every location there is a sheet added.
    for loc_key in str_dict:
        sheet = wb.add_worksheet(loc_key)
        form = wb.add_format()

        # The header is written on each sheets top.
        for col_index in range(len(header)):
            # Filter the argument from columns like feature(SIGNAL) or linage(PHYLUM).
            if re.findall(r'\(', header[col_index]):
                if re.findall(r'\((.*)\)', header[col_index]):
                    string = re.findall(r'\((.*)\)', header[col_index])[0]
                    if (string.lower() != "signal") or ("location" not in string.lower()):
                        string = "lineage"
            else:
                string = header[col_index]
            form.set_bold()
            sheet.write(0, col_index, string.upper(), form)

        # Calculate 3rd quantile of epitope density.
        if method.lower() != "uniprot":
            try:
                dens_string = list(list(zip(*list(str_dict[loc_key].values())))[1])
                densities = remove_nas(dens_string)
                q3_density = np.quantile(densities, .75)
            except IndexError:
                print(gene_id, "No Densities due to index error.")
        else:
            q3_density = None

        # Each line contains data of one protein that is indexed by ID.
        for gene_id in str_dict[loc_key]:
            form = wb.add_format()
            row_index = list(str_dict[loc_key].keys()).index(gene_id)

            # Color ID if gene density is above Q3 of all epitope densities in this location.
            if q3_density:
                try:
                    if float(str_dict[loc_key][gene_id][1]) > q3_density:
                        form.set_bg_color('#00FF00')
                except (TypeError, ValueError) as err:
                    print(gene_id, "No Density due to type or value error")
                    print(err)

            # Write ID as first column.
            sheet.write(row_index + 1, 0, gene_id, form)

            # Subsequently the columns are filled with data stored at gene ID.
            result_len = len(str_dict[loc_key][gene_id])
            link = 0
            for col_index in range(result_len):
                form = wb.add_format()
                if method.lower() == "uniprot":

                    # Protein names in italic.
                    if header[col_index + 1] == 'protein':
                        form.set_italic()
                    # Numbers in number format.
                    if header[col_index + 1] == ('length' or 'mass' or 'SIGNAL_SEQUENCE'):
                        form.set_num_format('0')
                    # Write Link at the end of each row in link format.
                    if header[col_index + 1] == 'link':
                        sheet.write_url(row_index + 1, col_index + 1, base + gene_id, string='UniProt')
                        link += 1
                else:
                    # Protein names in italic.
                    if col_index == 0:
                        form.set_italic()
                    # Numbers in number format.
                    if col_index in {1, 2}:
                        form.set_num_format('0.##0')
                    # Add border to epitopes
                    # Write link in link format and add border to epitope data.
                    if col_index == 4:
                        form.set_left()
                        sheet.write_url(row_index + 1, col_index + 1, base + gene_id, string='UniProt')
                        link += 1
                    # Write url to uniprot entry.
                    if col_index == 2 and result_len == 3:
                        sheet.write_url(row_index + 1, col_index + 2, base + gene_id, string='UniProt')
                    # Mark the epitopes inside signal.
                    if in_signal and col_index in in_signal[loc_key][gene_id]:
                        form.set_font_color('green')

                if isinstance(str_dict[loc_key][gene_id][col_index], int):
                    sheet.write(row_index + 1, col_index + 1 + link, int(str_dict[loc_key][gene_id][col_index]), form)
                elif isinstance(str_dict[loc_key][gene_id][col_index], float):
                    sheet.write(row_index + 1, col_index + 1 + link, float(str_dict[loc_key][gene_id][col_index]), form)
                elif str_dict[loc_key][gene_id][col_index].isnumeric():
                    sheet.write(row_index + 1, col_index + 1 + link, int(str_dict[loc_key][gene_id][col_index]), form)
                else:
                    sheet.write(row_index + 1, col_index + 1 + link, str_dict[loc_key][gene_id][col_index], form)

    wb.close()
    return 0


# Deletes redundant sequences from two dictionaries. In particular, the one printed to spreadsheets and the one given
# to the predition algorithms.
def del_redundant(data_dict, redundant_ids):
    popped = {}

    for key in list(data_dict):
        if key in redundant_ids:
            sequence = data_dict[key][-1]
            popped[key] = sequence
            data_dict.pop(key)
    return popped


# Set the parameters for predictions.
def pred_params(bool):
    if bool:
        allel = "HLA-A*02:01"
        length = "9"
        nmp = 3
        spt = 20
    else:
        allel = input("Allel: ")
        length = input("Length: ")
        print("Scores: highest -> lowest")
        spt = int(input("SYFPEITHI-Threshold [~36 -> 0]: "))
        nmp = float(input("NetMHCpan-Threshold [0.0 -> 100.0]: "))
    return allel, length, nmp, spt


# Print deleted sequences.
def print_popped(seq_dict):
    answer = input("Print popped sequences? [yes/no] ")
    if answer.lower() == "yes":
        deleted = 0
        for loc_key in seq_dict:
            print("\n" + loc_key + "\n")
            for gene_id in seq_dict[loc_key]:
                print(gene_id, "\t" + seq_dict[loc_key][gene_id] + "\n")
                deleted += 1
        print("Overall deleted:", deleted)
        return 0
    else:
        return 0


# Remove N/A entries from a list.
def remove_nas(str_list):
    return [element for element in str_list if element != 'N/A']

    return 0
