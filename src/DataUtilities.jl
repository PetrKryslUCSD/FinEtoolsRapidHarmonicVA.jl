using JSON
using DelimitedFiles
using HDF5
using FinEtools

const _sep = "-"

function getcleanname(name)
    cleanname = replace(replace(replace(name,'.'=>'_'),':'=>'_'),' '=>'_')    
    return cleanname
end

function with_extension(filename, ext)
    if match(Regex(".*\\." * ext * "\$"), filename) == nothing
        filename = filename * "." * ext
    end
    return filename
end

function retrieve_json(j)
    j = with_extension(j, "json")
    return open(j, "r") do file
        JSON.parse(file)
    end
end

function store_json(j, d)
    j = with_extension(j, "json")
	open(j, "w") do file
		JSON.print(file, d, 4)
	end
end

"""
    retrieve_matrix(fname)

Retrieve a matrix. 
"""
function retrieve_matrix(fname)
    fname = with_extension(fname, "h5")
    try
        return h5open(fname, "r") do file
            read(file, "matrix")
        end
    catch SystemError
        return nothing
    end
end

"""
    store_matrix(fname, matrix)

Store a matrix.
"""
function store_matrix(fname, matrix)
    fname = with_extension(fname, "h5")
    h5open(fname, "w") do file
        write(file, "matrix", matrix) 
    end
end
