


def to_fasta(df, outname, nuc_col_idx = 3, header_str = "1",
        other_data_idx = None):
    """
    write data as fasta format
    Returns: None
    """

    with open(outname, "w") as fout:
        hdr = ">" + header_str
        seq = "".join(df.iloc[:, nuc_col_idx].values)
        
        fout.write(hdr + "\n")
        fout.write(seq + "\n")
        
        if other_data_idx is not None:
            aux_dat = df.iloc[:, other_data_idx].values
            aux_dat = [str(x) for x in aux_dat]
            aux_dat = "".join(aux_dat)
            fout.write(aux_dat + "\n")
