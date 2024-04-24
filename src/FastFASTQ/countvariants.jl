
function countVariants(filepath::String,
                       ref::String; # reference sequence
                       seq_len::Integer=length(ref), 
                       start_pos::Integer=1,
                       min_qscore::Number=46,
                       min_avg_qscore::Number=60,
                       min_seq_len::Number=start_pos+seq_len+5,
                       maxMuts::Integer=2)
    # WARNING: THIS FUNCTION HAS NOT BEEN CHECKED FOR CORRECTNESS.
   
    # Count the number of variants in a Fastq file. 
    counts = zeros(Int, Int(seq_len/3), Int(seq_len/3), NUM_AA, NUM_AA) # 21 amino acids 
    ref = codeunits(ref) 
    muts = zeros(Int, maxMuts)
    pos = zeros(Int, maxMuts)

    fil = open(filepath, "r")
    i, N_qualReads, N_totHiMuts = 0, 0, 0
    while !eof(fil) 
        header = readline(fil) # read the header line
        seq = codeunits(readline(fil)) # read the sequence line
        plus = readline(fil) # read the "+" line 
        qual = codeunits(readline(fil)) # read the quality score line.
        if QC_read(qual, seq_len, start_pos, min_qscore, min_avg_qscore, min_seq_len)
            N_qualReads += 1 
            updateCounts!(counts, seq, ref, muts, pos, seq_len, start_pos) && (N_totHiMuts += 1)
        end
        i+=1
    end
    close(fil)
    N_reads = i 
    return counts, N_reads, N_qualReads, N_totHiMuts
end

function QC_read(qual::Base.CodeUnits, seq_len, start_pos, min_qscore, min_avg_qscore, min_seq_len)
    # return false if any of the following are satisfied:
    # 1. Any qscore is less than min_qscore.
    # 2. The average qscore is less then min_avg_qscore.
    # 3. The length of the read is less than min_seq_len
    
    if length(qual) < min_seq_len
        return false
    end
    sum_q = zero(Int)
    @inbounds for i in 1:seq_len
        q = Int(qual[i + start_pos - 1])
        if q < min_qscore
            return false
        end
        sum_q += q
    end
    avg_q = sum_q / seq_len
    if avg_q < min_avg_qscore
        return false
    else
        return true
    end
end


function updateCounts!(counts, seq::Base.CodeUnits, ref::Base.CodeUnits, muts, pos,
                       seq_len, start_pos; maxMuts=2)
    # return the mutant indicies
    numMuts = 0
    seq_len_AA = Int(seq_len/3) 
    @inbounds for i in 1:seq_len_AA
        s1 = seq[start_pos + 3i-3] 
        s2 = seq[start_pos + 3i-2]
        s3 = seq[start_pos + 3i-1]
        r1 = ref[3i-2]
        r2 = ref[3i-1]
        r3 = ref[3i  ]
        AA_seq = TABLE[codon2Byte(s1,s2,s3)]
        AA_ref = TABLE[codon2Byte(r1,r2,r3)]
        # check if mutated
        if (AA_seq != AA_ref) && numMuts <= maxMuts  # if AA mutated
            numMuts += 1
            muts[numMuts] = Int(AA_seq)
            pos[numMuts] = i
        end
    end

    higherMut = false
    if numMuts<1
        nothing
    elseif numMuts==1
        i = pos[1] 
        AA_i = muts[1]
        counts[i,i, AA_i, AA_i] +=1 
    elseif numMuts==2
        i,j = pos 
        AA_i, AA_j = muts
        counts[i,j, AA_i, AA_j] +=1 
    else
        higherMut = true
    end
    return higherMut
end





