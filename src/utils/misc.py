def _infer_bin_prefix(bin_name):
    prefix = ''
    for i in range(len(bin_name)):
        if bin_name[i].isalpha():
            prefix += bin_name[i]
        else:
            break
    return prefix


def correct_bin_names(bin_names):
    for i in range(len(bin_names)):
        bin_names[i] = bin_names[i].replace(':', '_', bin_names[i].count(':')-1)
    return bin_names