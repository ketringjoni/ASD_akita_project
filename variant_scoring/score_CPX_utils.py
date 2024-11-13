




def get_sequence_CPX(CHR, POS, END, CHR2, CPX_TYPE, CPX_INTERVALS):
    
    
    # Get reference and alternate sequence from REF and ALT allele using reference genome 


    # Get reference sequence

    REF_len = END - POS

    assert(REF_len < MB)

    REF_half_left = math.ceil((MB - REF_len)/2) # if the REF allele is odd, shift right
    REF_half_right = math.floor((MB - REF_len)/2)

    # Annotate whether variant is close to beginning or end of chromosome
    var_position = get_variant_position(CHR, POS, REF_len, REF_half_left, REF_half_right)

    # Get last coordinate of chromosome
    chrom_max = int(hg38_lengths[hg38_lengths.CHROM == CHR[3:]]['chrom_max']) 

    # Get centromere coordinate
    centromere = int(centromere_coords[centromere_coords.CHROM == CHR]['centromere'])


    # Get start and end of reference sequence

    if var_position == "chrom_mid":
        REF_start = POS - REF_half_left
        REF_stop = REF_start + MB 

    elif var_position == "chrom_start": 
        REF_start = 0
        REF_stop = MB

    elif var_position == "chrom_centro_left": 
        REF_start = centromere - MB
        REF_stop = centromere

    elif var_position == "chrom_centro_right": 
        REF_start = centromere
        REF_stop = centromere + MB

    elif var_position == "chrom_end": 
        REF_start = chrom_max - MB
        REF_stop = chrom_max


    # Get reference sequence
    REF_seq = fasta_open.fetch(CHR, REF_start, REF_stop).upper()


    # Error if Ns are more than 5% of sequence
    if Counter(REF_seq)['N']/MB*100 > 5:
        raise ValueError('N composition greater than 5%')



    # Get alternate sequence


    # get sequence for inversion

    # get sequence to invert
    INV_coordinates = CPX_INTERVALS.split('INV_' + CHR2.split('chr')[1] + ':')[1].split(',')[0]
    INV_start = int(INV_coordinates.split('-')[0])
    INV_end = int(INV_coordinates.split('-')[1])

    inv_revcomp = fasta_open.fetch(CHR2, INV_start, INV_end).upper()
    inv_seq = str(Seq(inv_revcomp).reverse_complement())



    # get sequences left and right of inversion


    # Get variant sequence before inversion
    if CPX_TYPE in ['delINVdel', 'INVdup', 'delINVdup']:
        ALT_left = ''

    if CPX_TYPE in ['dupINVdup', 'dupINV', 'dupINVdel']:
        DUP1_coordinates = CPX_INTERVALS.split('DUP_' + CHR2.split('chr')[1] + ':')[1].split(',')[0]
        DUP1_start = int(DUP1_coordinates.split('-')[0])
        DUP1_end = int(DUP1_coordinates.split('-')[1])
        ALT_left = fasta_open.fetch(CHR2, DUP1_start, DUP1_end).upper()


    # Get variant sequence after inversion
    if CPX_TYPE in ['delINVdel', 'dupINV', 'dupINVdel']:
        ALT_right = ''

    if CPX_TYPE in ['dupINVdup', 'INVdup', 'delINVdup']:

        if CPX_TYPE != 'dupINVdup':
            num = 1
        # if there are 2 duplications, get coordinates for the seconds one
        else:
            num = 2

        DUP2_coordinates = CPX_INTERVALS.split('DUP_' + CHR2.split('chr')[1] + ':')[num].split(',')[0]
        DUP2_start = int(DUP2_coordinates.split('-')[0])
        DUP2_end = int(DUP2_coordinates.split('-')[1])

        ALT_right = fasta_open.fetch(CHR2, DUP2_start, DUP2_end).upper()


    # Assemble ALT variant sequence
    ALT_var_seq = ALT_left + inv_seq + ALT_right

    assert(len(ALT_var_seq) < MB)

    # Get ALT left and right sequences from ref genome
    ALT_half_left = math.ceil((MB - len(ALT_var_seq))/2) # if the ALT allele is odd, shift right
    ALT_half_right = math.floor((MB - len(ALT_var_seq))/2)

    # Get start and end of reference sequence

    if var_position == "chrom_mid":
        ALT_leftREF = fasta_open.fetch(CHR, POS - ALT_half_left, POS).upper()
        ALT_rightREF = fasta_open.fetch(CHR, END, END + ALT_half_right).upper()

    elif var_position == "chrom_start": 
        ALT_leftREF = fasta_open.fetch(CHR, 0, POS).upper()
        ALT_half_right = MB - POS - len(ALT_var_seq)
        ALT_rightREF = fasta_open.fetch(CHR, END, END + ALT_half_right).upper()

    elif var_position == "chrom_centro_left": 
        ALT_rightREF = fasta_open.fetch(CHR, END, centromere).upper()
        ALT_half_left = MB - (centromere - END) - len(ALT_var_seq)
        ALT_leftREF = fasta_open.fetch(CHR, POS - ALT_half_left, POS).upper()

    elif var_position == "chrom_centro_right": 
        ALT_leftREF = fasta_open.fetch(CHR, centromere, POS).upper()
        ALT_half_right = MB - (POS - centromere) - len(ALT_var_seq)
        ALT_rightREF = fasta_open.fetch(CHR, END, END + ALT_half_right).upper()

    elif var_position == "chrom_end": 
        ALT_rightREF = fasta_open.fetch(CHR, END, chrom_max).upper()
        ALT_half_left = MB - (chrom_max - END) - len(ALT_var_seq)
        ALT_leftREF = fasta_open.fetch(CHR, POS - ALT_half_left, POS).upper()



    # Assemble ALT sequence
    ALT_seq = ALT_leftREF + ALT_var_seq + ALT_rightREF


    # make sure ALT_seq is correct length
    assert(len(ALT_seq) == MB)



    return REF_seq, ALT_seq






