import zipfile
import glob


def extract_files(basezip, outzip, string, indicator='effpre', verbose=True):
    with zipfile.ZipFile(basezip, 'r') as base_archive:
        filelist = base_archive.namelist()
        files = [file for file in filelist if string in file]
        files = [file for file in files if indicator in file]

        with zipfile.ZipFile(outzip, 'w', compression=zipfile.ZIP_DEFLATED) as out_archive:
            for file in files:
                file_contents = base_archive.read(file)
                out_archive.writestr(file, file_contents)
                if verbose:
                    print(f"saved {file}")


if __name__ == "__main__":
    basezip = '/Users/alison/Downloads/w@h effective precip.zip'                    
    extract_files(basezip, "/Users/alison/Documents/RAPID/data_input/w@h/bs.zip", "baseline")
    extract_files(basezip, "/Users/alison/Documents/RAPID/data_input/w@h/nf.zip", "near_future")
    extract_files(basezip, "/Users/alison/Documents/RAPID/data_input/w@h/ff.zip", "far_future")