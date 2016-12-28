using DataFrames
export load_narrowpeak, metadata

metadata = readtable(joinpath(Pkg.dir("ChipSeqUtil"), "src", "metadata.csv"))

function load_narrowpeak(stream, contigs, index; binSize=1000, loadp=true, mlogt=false, verbose=0)
    numBins = ceil(Int64, sum(contigs.sizes) / binSize)
    chrOffsets = Dict{String,Int64}()
    for i in 1:contigs.count
        chrOffsets[contigs.names[i]] = contigs.offsets[i]
    end

    # mark all bins that are touched with 1
    binValues = falses(numBins)
    pValues = zeros(numBins)
    count = 0 
    for line in eachline(stream)
        parts = split(line, '\t')
        if haskey(chrOffsets, parts[1])
            count += 1
            if count < verbose
                println(parts)
            end
            startPos = ceil(Int64, (chrOffsets[parts[1]]+parse(Int64, parts[2]))/binSize)
            endPos = ceil(Int64, (chrOffsets[parts[1]]+parse(Int64, parts[3]))/binSize)
            for i in startPos:endPos
                if i!=0
                    binValues[i] = true
                    if loadp
                        if mlogt
                            pval = -1*log(10, parse(Float64, parts[index]))
                        else
                            pval = parse(Float64, parts[index])
                        end
                        if pValues[i] < pval
                            pValues[i] = pval
                        end
                    end
                end
            end
        end
    end
    close(stream)
    binValues, pValues
end
