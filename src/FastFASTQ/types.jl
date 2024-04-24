

mutable struct SequenceRead
    endIndex::Int
    seqBuffer::Vector{UInt8}
    qualBuffer::Vector{UInt8}
end

mutable struct BlockReader{IOT <: IO}
    io::IOT
    byteVec::Vector{UInt8}
    nb::Int # Number of bytes read by readbytes
    startFirstFullRead::Int
    endLastFullRead::Int
    currentIndex::Int
end



