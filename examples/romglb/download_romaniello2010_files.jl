import Downloads
import ZipFile

function download_romaniello2010_files(;
    url="https://github.com/PALEOtoolkit/PALEOocean.jl/releases/download/v0.4.6/romaniello2010_transport.zip",
    output_folder=@__DIR__,
)

    # download zip archive to temporary location
    zarchive_path = Downloads.download(url)

    # unpack into a new subfolder
    zarchive = ZipFile.Reader(zarchive_path)
    try
        for f in zarchive.files
            full_file_path = joinpath(output_folder, f.name)
            if (endswith(f.name,"/") || endswith(f.name,"\\"))
                mkdir(full_file_path)
            else
                write(full_file_path, read(f))
            end
        end
    finally
        close(zarchive)
    end

    return nothing
end