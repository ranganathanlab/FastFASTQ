function updateCounts!(counts, seq::Base.CodeUnits, ref::Base.CodeUnits, AATable::Array,
                       seq_len, start_pos; maxMuts=3)
    # return the mutant indicies
    muts = ()
    mut_tally = 0
    seq_len_AA = Int(seq_len/3) 
    for i in 1:seq_len_AA
        s1 = Int(seq[start_pos + 3i-3]) -  ORDSHIFT 
        s2 = Int(seq[start_pos + 3i-2]) -  ORDSHIFT
        s3 = Int(seq[start_pos + 3i-1]) -  ORDSHIFT
        r1 = Int(ref[3i-2]) -  ORDSHIFT
        r2 = Int(ref[3i-1]) -  ORDSHIFT
        r3 = Int(ref[3i  ]) -  ORDSHIFT
        # check if mutated
        if s1 != r1 || s2 != r2 || s3 != r3 # if DNA mutated
            AA_seq = AATable[s1,s2,s3]
            AA_ref = AATable[r1,r2,r3]
            if AA_seq != AA_ref # if AA mutated
                muts = (muts..., (i, AA_seq)) 
            end
        end
    end

    higherMut = false
    numMuts = length(muts)
    if numMuts<1
        nothing
    elseif numMuts==1
        i, AA_i = muts[1]
        counts[i,i, AA_i, AA_i] +=1 
    elseif numMuts==2
        i, AA_i = muts[1]
        j, AA_j = muts[2]
        counts[i,j, AA_i, AA_j] +=1 
    else
        higherMut = true
    end

    return higherMut
end

function buildAATable()
    # build an 3d Int array whose index i,j,k specify a codon
    # and the value at i,j,k is the amino acid integer label.
    AATable = zeros(Int, 20, 20, 20)
    for codon in keys(CODON2AA)
        aa = CODON2AA[codon]
        c1,c2,c3 = Int.(codeunits(codon)) .- ORDSHIFT
        AATable[c1,c2,c3] = AA2INT[aa]
    end
    return AATable
end

