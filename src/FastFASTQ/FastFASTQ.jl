module FastFASTQ

# using

push!(LOAD_PATH,"./")
include("types.jl")
include("countvariants.jl")
#include("new_code.jl")


##### Export Functions ##################################
export 

    # Constants ----------------------------------------
    CODON2AA,
    AA2INT,
    TABLE,
    
    # Types
    BlockReader,
    SequenceRead,

    # Functions ---------------------------------------
    codon2Byte,

    # from countvariants.jl
    countVariants,
    QC_read,
    updateCounts!,

    # from new_code.jl
    countVariants2,
    countVariants3,
    getNextRead!,
    updateBlockInfo!,
    get_startFirstFullRead,
    get_endLastFullRead,
    getNextBlock!,
    buildBlockReader




##### Define Constants ##################################

const NUM_AA::Int = 21 # Number of Amino Acids including Stop codon.

const CODON2AA::Dict{String, String} = Dict( 
    "ACC" => "T", "ATG" => "M", "ACA" => "T",
    "ACG" => "T", "ATC" => "I", "AAC" => "N",
    "ATA" => "I", "AGG" => "R", "CCT" => "P",
    "CTC" => "L", "AGC" => "S", "AAG" => "K",
    "AGA" => "R", "CAT" => "H", "AAT" => "N",
    "ATT" => "I", "CTG" => "L", "CTA" => "L",
    "ACT" => "T", "CAC" => "H", "AAA" => "K",
    "CCG" => "P", "AGT" => "S", "CCA" => "P",
    "CAA" => "Q", "CCC" => "P", "TAT" => "Y",
    "GGT" => "G", "TGT" => "C", "CGA" => "R",
    "CAG" => "Q", "CGC" => "R", "GAT" => "D",
    "CGG" => "R", "CTT" => "L", "TGC" => "C",
    "GGG" => "G", "TAG" => "*", "GGA" => "G",
    "TAA" => "*", "GGC" => "G", "TAC" => "Y",
    "GAG" => "E", "TCG" => "S", "TTT" => "F",
    "GAC" => "D", "CGT" => "R", "GAA" => "E",
    "TCA" => "S", "GCA" => "A", "GTA" => "V",
    "GCC" => "A", "GTC" => "V", "GCG" => "A",
    "GTG" => "V", "TTC" => "F", "GTT" => "V",
    "GCT" => "A", "TTA" => "L", "TGA" => "*",
    "TTG" => "L", "TCC" => "S", "TGG" => "W",
    "TCT" => "S")

const AA2INT::Dict{String, Int} = Dict(
     "A"=>1,"C"=>2, "D"=>3, "E"=>4, "F"=>5, "G"=>6, "H"=>7,
     "I"=>8, "K"=>9, "L"=>10, "M"=>11, "N"=>12, "P"=>13, "Q"=>14,
     "R"=>15, "S"=>16, "T"=>17, "V"=>18, "W"=>19, "Y"=>20, "*"=>21)

function codon2Byte(c1::UInt8, c2::UInt8, c3::UInt8)
    # returns a unique byte for each codon
    d1 = c1 >> 1 << 6 >> 6
    d2 = c2 >> 1 << 6 >> 4
    d3 = c3 >> 1 << 6 >> 2
    return d1 + d2 + d3 + one(UInt8) # the one is added for indexing the TABLE
end

const TABLE::NTuple{64, UInt8} = let  
    table = fill(0xff, 2^6) # 64 long
    for codon in keys(CODON2AA)
        aa = CODON2AA[codon]
        aa_int = AA2INT[aa]
        c1,c2,c3 = codeunits(codon)
        byte = codon2Byte(c1,c2,c3)
        table[byte] = UInt8(aa_int)
    end
    Tuple(table)
end


end # end of module
































#function main(filepath::String,
#              ref::String;
#              seq_len=length(ref), 
#              start_pos=1,
#              min_qscore=46,
#              min_avg_qscore=60,
#              min_seq_len=start_pos+seq_len+5,
#              maxMuts = 2)
#   
#    # process Fastq file.
#    
#    AATable = buildAATable() # create AA table
#    counts = zeros(Int, Int(seq_len/3), Int(seq_len/3), NUM_AA, NUM_AA) # 21 amino acids 
#    ref = codeunits(ref) 
#    muts = zeros(Int, maxMuts)
#    pos = zeros(Int, maxMuts)
#
#    fil = open(filepath, "r")
#    i, t, totalHigherMuts = 0, 0, 0
#    while !eof(fil)
#        header = readline(fil)
#        seq = codeunits(readline(fil))
#        plus = readline(fil)
#        qual = codeunits(readline(fil))
#        if QC_read(qual, seq_len, start_pos, min_qscore, min_avg_qscore, min_seq_len)
#            t+=1
#            updateCounts2!(counts, seq, ref, muts, pos, seq_len, start_pos) && (totalHigherMuts += 1)
#        end
#        i+=1
#    end
#    close(fil)
#    return counts, t, i, totalHigherMuts
#end

#function QC_read(qual::Base.CodeUnits, seq_len, start_pos, min_qscore, min_avg_qscore, min_seq_len)
#    # return false if any of the following are satisfied:
#    # 1. Any qscore is less than min_qscore.
#    # 2. The average qscore is less then min_avg_qscore.
#    # 3. The length of the read is less than min_seq_len
#    
#    if length(qual) < min_seq_len
#        return false
#    end
#    sum_q = zero(Int)
#    @inbounds for i in 1:seq_len
#        q = Int(qual[i + start_pos - 1])
#        if q < min_qscore
#            return false
#        end
#        sum_q += q
#    end
#    avg_q = sum_q / seq_len
#    if avg_q < min_avg_qscore
#        return false
#    else
#        return true
#    end
#end

#function updateCounts!(counts, seq::Base.CodeUnits, ref::Base.CodeUnits, AATable::Array,
#                       seq_len, start_pos; maxMuts=3)
#    # return the mutant indicies
#    muts = ()
#    mut_tally = 0
#    seq_len_AA = Int(seq_len/3) 
#    for i in 1:seq_len_AA
#        s1 = Int(seq[start_pos + 3i-3]) -  ORDSHIFT 
#        s2 = Int(seq[start_pos + 3i-2]) -  ORDSHIFT
#        s3 = Int(seq[start_pos + 3i-1]) -  ORDSHIFT
#        r1 = Int(ref[3i-2]) -  ORDSHIFT
#        r2 = Int(ref[3i-1]) -  ORDSHIFT
#        r3 = Int(ref[3i  ]) -  ORDSHIFT
#        # check if mutated
#        if s1 != r1 || s2 != r2 || s3 != r3 # if DNA mutated
#            AA_seq = AATable[s1,s2,s3]
#            AA_ref = AATable[r1,r2,r3]
#            if AA_seq != AA_ref # if AA mutated
#                muts = (muts..., (i, AA_seq)) 
#            end
#        end
#    end
#
#    higherMut = false
#    numMuts = length(muts)
#    if numMuts<1
#        nothing
#    elseif numMuts==1
#        i, AA_i = muts[1]
#        counts[i,i, AA_i, AA_i] +=1 
#    elseif numMuts==2
#        i, AA_i = muts[1]
#        j, AA_j = muts[2]
#        counts[i,j, AA_i, AA_j] +=1 
#    else
#        higherMut = true
#    end
#
#    return higherMut
#end
#
#function updateCounts2!(counts, seq::Base.CodeUnits, ref::Base.CodeUnits, muts, pos,
#                       seq_len, start_pos; maxMuts=2)
#    # return the mutant indicies
#    numMuts = 0
#    seq_len_AA = Int(seq_len/3) 
#    @inbounds for i in 1:seq_len_AA
#        s1 = seq[start_pos + 3i-3] 
#        s2 = seq[start_pos + 3i-2]
#        s3 = seq[start_pos + 3i-1]
#        r1 = ref[3i-2]
#        r2 = ref[3i-1]
#        r3 = ref[3i  ]
#        AA_seq = TABLE[codon2Byte(s1,s2,s3)]
#        AA_ref = TABLE[codon2Byte(r1,r2,r3)]
#        # check if mutated
#        if (AA_seq != AA_ref) && numMuts <= maxMuts  # if AA mutated
#            numMuts += 1
#            muts[numMuts] = Int(AA_seq)
#            pos[numMuts] = i
#        end
#    end
#
#    higherMut = false
#    if numMuts<1
#        nothing
#    elseif numMuts==1
#        i = pos[1] 
#        AA_i = muts[1]
#        counts[i,i, AA_i, AA_i] +=1 
#    elseif numMuts==2
#        i,j = pos 
#        AA_i, AA_j = muts
#        counts[i,j, AA_i, AA_j] +=1 
#    else
#        higherMut = true
#    end
#    return higherMut
#end

######################################################





#function buildAATable()
#    # build an 3d Int array whose index i,j,k specify a codon
#    # and the value at i,j,k is the amino acid integer label.
#    AATable = zeros(Int, 20, 20, 20)
#    for codon in keys(CODON2AA)
#        aa = CODON2AA[codon]
#        c1,c2,c3 = Int.(codeunits(codon)) .- ORDSHIFT
#        AATable[c1,c2,c3] = AA2INT[aa]
#    end
#    return AATable
#end


