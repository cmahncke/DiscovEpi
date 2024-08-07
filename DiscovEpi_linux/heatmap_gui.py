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

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pylab as plt
import time
import addtional_guis as add


def produce_heatmap(prot_res, epis, dir, topnr, cutoff, progress, worker):
    try:
        start = time.time()
        plt.style.use("seaborn-v0_8")
        prot = prot_res[0]
        namestem = prot_res[5]
        epis = epis[0]
        cutoff = int(cutoff) if cutoff else 1000000
        topnr = int(topnr) if topnr else 1000000

        for loc_key in prot:

            if worker.isInterruptionRequested():
                return -2

            unp = prot[loc_key]
            nmp = epis[loc_key]
            nmp = dict(sorted(nmp.items(), key=lambda item: item[1][1], reverse = True))
            key_list = list(nmp.keys())[:topnr] if topnr < len(nmp.keys()) else list(nmp.keys())

            top10max = 0
            max_value = 0
            hm_temp_scores_data = np.array([[0.0]])
            max_bind = 0
            index = 0

            for gene in nmp:

                if worker.isInterruptionRequested():
                    return -2

                index += 1
                if index > topnr: break

                seq = unp[gene][0]
                seq = seq[:cutoff] if len(seq) > cutoff else seq
                length = len(seq)

                pos_scores_lst = [0.5] * length

                for epitope in nmp[gene][4:]:

                    if worker.isInterruptionRequested():
                        return -2

                    if isinstance(epitope, str):
                        pos = int(epitope.split(";")[0])
                        if pos+8 > cutoff: continue
                        bind_score = float(epitope.split(";")[2])
                        for i in range(pos, pos + 8):
                            pos_scores_lst[i - 1] += bind_score

                rest = np.shape(hm_temp_scores_data)[1] - len(pos_scores_lst)
                if rest < 0:
                    appendix = [[0.0] * abs(rest)] * np.shape(hm_temp_scores_data)[0]
                    hm_temp_scores_data = np.append(hm_temp_scores_data, list(appendix), axis=1)
                    hm_temp_scores_data = np.append(hm_temp_scores_data, [pos_scores_lst], axis=0)

                elif rest > 0:
                    appendix = [0.0] * rest
                    pos_scores_lst += appendix
                    hm_temp_scores_data = np.append(hm_temp_scores_data, [pos_scores_lst], axis=0)

                else:
                    hm_temp_scores_data = np.append(hm_temp_scores_data, [pos_scores_lst], axis=0)

                if max_bind < max(pos_scores_lst):
                    max_bind = max(pos_scores_lst)

                progress.emit((list(key_list).index(gene)+1) / len(key_list) * 100)

            if worker.isInterruptionRequested():
                return -2

            hm_temp_scores_data = np.delete(hm_temp_scores_data, 0, axis=0)
            shape = np.shape(hm_temp_scores_data)
            hm_score_data = np.array([[0.0]*shape[1]]*shape[0])
            for ix1, protein in enumerate(hm_temp_scores_data):
                for ix2, score in enumerate(protein):
                    hm_score_data[ix1][ix2] = score
                    if score > max_value:
                        max_value = score
                    if ix1 < 100 and score > top10max:
                        top10max = score

            hm_score_data = hm_score_data
            hm_df = pd.DataFrame(hm_score_data, index=key_list)
            plt.figure(figsize=(30, 20))
            sns.set(font_scale=3)

            if worker.isInterruptionRequested():
                return -2

            heat_map = sns.heatmap(hm_df, vmin=0, vmax=max_value, xticklabels = 100, cmap=plt.colormaps["Greys"], cbar=False)

            filename = f"hmp_{namestem}_{loc_key}.png"
            plt.title("Epitope map: " + loc_key, fontsize = 40)
            plt.xlabel("Amino acid sequence [absolute position of amino acid]", fontsize = 40)
            path = dir + "\\" + filename
            print(path)
            plt.savefig(path)

        end = time.time()
        print("Time elapsed in heatmap: ", add.calc_time(end - start))

        return 0
    except:
        print("Something wrong with the heatmap.")
        return -1