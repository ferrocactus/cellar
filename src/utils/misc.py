def _infer_bin_prefix(bin_name):
    prefix = ''
    for i in range(len(bin_name)):
        if bin_name[i].isalpha():
            prefix += bin_name[i]
        else:
            break
    return prefix
