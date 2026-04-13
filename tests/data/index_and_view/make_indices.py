import pickle


def make_index():
    node_info = {}

    node_info["s1"] = ["REF#0#CONTIG1", 0, 10]
    node_info["s2"] = ["REF#0#CONTIG1", 10, 30]
    node_info["s3"] = ["REF#0#CONTIG1", 30, 35]
    node_info["s4"] = ["REF#0#CONTIG1", 35, 60]
    node_info["s14"] = ["FOO#1#ASSM1", 10, 30]
    node_info["s28"] = ["BAR#1#ASSM1", 240, 250]
    node_info["s6"] = ["REF#0#CONTIG1", 90, 100]
    node_info["s15"] = ["FOO#1#ASSM1", 60, 85]
    node_info["s16"] = ["FOO#1#ASSM1", 85, 90]
    node_info["s22"] = ["FOO#2#ASSM2", 80, 85]
    node_info["s27"] = ["BAR#1#ASSM1", 200, 230]
    node_info["s32"] = ["BAZ#1#ASSM1", 105, 125]
    node_info["s8"] = ["REF#0#CONTIG1", 110, 125]
    node_info["s20"] = ["FOO#2#ASSM1", 175, 205]
    node_info["s18"] = ["FOO#1#ASSM1", 175, 200]
    node_info["s12"] = ["REF#0#CONTIG1", 180, 200]
    node_info["s29"] = ["BAR#2#ASSM1", 10, 30]
    node_info["s24"] = ["BAR#1#ASSM1", 10, 40]
    node_info["s21"] = ["FOO#2#ASSM1", 205, 215]
    node_info["s13"] = ["REF#0#CONTIG1", 200, 210]

    # re-formating the node_info dict
    for id in node_info.keys():
        sn, so, se = node_info[id]
        node_info[id] = (id, sn, so, se)

    stable_offsets = {}

    stable_offsets["s1"] = [0, 97, 1545, 1779]
    stable_offsets["s2"] = [0]
    stable_offsets["s3"] = [0, 235, 1545, 1779]
    stable_offsets["s4"] = [0, 97, 235]
    stable_offsets["s14"] = [97, 1545, 1779]
    stable_offsets["s28"] = [235]
    stable_offsets["s6"] = [235, 409, 2212]
    stable_offsets["s15"] = [409, 2212]
    stable_offsets["s16"] = [409, 2212]
    stable_offsets["s22"] = [409, 2212]
    stable_offsets["s27"] = [552, 655, 823, 1926, 2100]
    stable_offsets["s32"] = [655, 1926]
    stable_offsets["s8"] = [655, 1926]
    stable_offsets["s20"] = [925, 1422]
    stable_offsets["s18"] = [1020, 1115, 1422]
    stable_offsets["s12"] = [1115, 1683]
    stable_offsets["s29"] = [1236]
    stable_offsets["s24"] = [1329]
    stable_offsets["s21"] = [1422]
    stable_offsets["s13"] = [1683]

    unstable_offsets = {}

    unstable_offsets["s1"] = [0, 95, 1201, 1379]
    unstable_offsets["s2"] = [0]
    unstable_offsets["s3"] = [0, 186, 1201, 1379]
    unstable_offsets["s4"] = [0, 95, 186]
    unstable_offsets["s14"] = [95, 1201, 1379]
    unstable_offsets["s28"] = [186]
    unstable_offsets["s6"] = [186, 312, 1698]
    unstable_offsets["s15"] = [312, 1698]
    unstable_offsets["s16"] = [312, 1698]
    unstable_offsets["s22"] = [312, 1698]
    unstable_offsets["s27"] = [413, 500, 617, 1479, 1602]
    unstable_offsets["s32"] = [500, 1479]
    unstable_offsets["s8"] = [500, 1479]
    unstable_offsets["s20"] = [703, 1106]
    unstable_offsets["s18"] = [782, 861, 1106]
    unstable_offsets["s12"] = [861, 1292]
    unstable_offsets["s29"] = [948]
    unstable_offsets["s24"] = [1027]
    unstable_offsets["s21"] = [1106]
    unstable_offsets["s13"] = [1292]

    stable_dict = {}
    for id in node_info.keys():
        stable_dict[node_info[id]] = stable_offsets[id]
    stable_dict["ref_contig"] = ["REF#0#CONTIG1"]

    with open("view-index-stable.gvi", "wb") as handle:
        pickle.dump(stable_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    unstable_dict = {}
    for id in node_info.keys():
        unstable_dict[node_info[id]] = unstable_offsets[id]
    unstable_dict["ref_contig"] = ["REF#0#CONTIG1"]

    with open("view-index-unstable.gvi", "wb") as handle:
        pickle.dump(unstable_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    make_index()
