

function updateBlockInfo!(block::BlockReader)
    # sets the Block index positions.
    block.startFirstFullRead = get_startFirstFullRead(block)
    block.endLastFullRead = get_endLastFullRead(block)
    block.currentIndex = 1
    return nothing
end


function get_startFirstFullRead(block::BlockReader; startReadSymbol::UInt8=UInt8('@'))
    # get the index of the start of the first full
    # read in the block
    for i in 1:block.nb 
        byte = block.byteVec[i]
        if byte == startReadSymbol
            return i
        end
    end
    println("No @ was found in this block. Must be the end")
    return -1 # indicates no read start found
end

function get_endLastFullRead(block::BlockReader; startReadSymbol::UInt8=UInt8('@'))
    # TEST EDGE CASES FOR THIS FUNCTION
    # get the index of the end of the last full
    # read of the block.
    nl = UInt8('\n')
    N_nl = 0 # Number of newlines counted
    nb = block.nb
    i = nb
    byte = block.byteVec[i]
    while byte != startReadSymbol && i > 1
        byte == nl && (N_nl += 1)
        i -= 1
        byte = block.byteVec[i]
    end
    
    if N_nl < 4
        return i
    elseif N_nl == 4
        return nb
    else
        error("The end of FASTQ not formatted properly!")
    end
end

function getNextBlock!(block::BlockReader)
    eof(block.io) && error("Trying to get Next block when EOF")
    nb = readbytes!(block.io, block.byteVec)
    nb < length(block.byteVec) && println("End of File")
    block.nb = nb
    updateBlockInfo!(block)
    return nothing
end


function buildBlockReader(filepath::String, blockSize::Int)
    @assert blockSize > 10_000
    return BlockReader(open(filepath, "r"), zeros(UInt8, blockSize),
                       blockSize, 1, blockSize, 1)
end

function countVariants3(filepath::String; blockSize=304_820)
    
    block = buildBlockReader(filepath, blockSize)
    while !eof(block.io)
        getNextBlock!(block)
       # @show block.nb
       # @show block.startFirstFullRead
       # @show block.endLastFullRead
       # println()
    end
    close(block.io)
    return nothing
end




function getNextRead!(block::BlockReader, startReadSymbol::UInt8=UInt8('@'))
    # extra is a vector that goes from the start
    # of the interupted previous read to the interupt
    # site

    startReadSymbol=UInt8('@')
    nl=UInt8('@')
    byteVec = block.byteVec
    current = block.currentIndex
    nls = zeros(Int, 4)
    N_nl = 0
    i = current
    c = byteVec[i]
    while c != startReadSymbol
         if c == nl
             N_nl+=1
             nls[N_nl] = i
         end
         i+=1
         c=
    end
    

end




function countVariants2(filepath::String,
                       ref::String;
                       seq_len=length(ref), 
                       start_pos=1,
                       min_qscore=46,
                       min_avg_qscore=60,
                       min_seq_len=start_pos+seq_len+5,
                       maxMuts = 2)
   
    # Count the number of variants in a Fastq file.
    counts = zeros(Int, Int(seq_len/3), Int(seq_len/3), NUM_AA, NUM_AA) # 21 amino acids 
    ref = codeunits(ref) 
    muts = zeros(Int, maxMuts)
    pos = zeros(Int, maxMuts)
    seqread = SequenceRead(1, zeros(UInt8, 2seq_len), zeros(UInt8,2seq_len))

    fil = open(filepath, "r")
    i, N_qualReads, N_totHiMuts = 0, 0, 0
    while !eof(fil)
        getNextRead!(seqread, fil)
        #goodRead = QC_read(qual, seq_len, start_pos, min_qscore, min_avg_qscore, min_seq_len)
        #if goodRead
            #N_qualReads += 1
            #updateCounts!(counts, seq, ref, muts, pos, seq_len, start_pos) && (N_totHiMuts += 1)
        #end
        i+=1
    end
    close(fil)
    N_reads = i 
    return counts, N_reads, N_qualReads, N_totHiMuts
end

function getNextRead!(seqread::SequenceRead, fastqIO::IO)
    # takes an io stream from a FASTQ formated file and 
    # reads out one FASTQ read, writing the sequence and
    # quality scores to buffers in the read object. 
    # This function can be called iteratively on a single
    # IO steam to read the whole file.
    newLineCount = 0
    writeSeq = writeQual = false
    i = seq_i = qual_i = 0
    while newLineCount < 4 # FASTQ format has 4 lines per read
        i += 1
        c = read(fastqIO, UInt8)

        if writeSeq && c!=0x0a
            seq_i += 1
            seqread.seqBuffer[seq_i] = c
        elseif writeQual && c!=0x0a
            qual_i += 1
            seqread.qualBuffer[qual_i] = c
        end

        if c == 0x0a # new line character
            if newLineCount==0 # first '\n'
                newLineCount += 1
                write_seq = true
                write_qual = false
            elseif newLineCount==1 # second '\n'
                newLineCount += 1
                seqread.endIndex = seq_i # record the end index of the seq.
                read(fastqIO, UInt8) # read the plus character
                read(fastqIO, UInt8) # read the third '\n'
                newLineCount += 1 
                writeSeq = false
                writeQual = true
            elseif newLineCount==3 # third '\n'
                newLineCount += 1
                writeSeq = false
                writeQual = false
            end
        end
    end
    return nothing
end
    

