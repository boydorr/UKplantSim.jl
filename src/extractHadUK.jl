using PyCall
py"""
# Import required python modules
import ftplib
import os

def retrieve_HadUK(param, from_year, to_year, directory, user = "charris009", passwd = "HadUK3019"):
    os.chdir(directory)
    # login to FTP
    f=ftplib.FTP("ftp.ceda.ac.uk", user, passwd)

    # loop through years
    for year in range(from_year,to_year+1):
      # change the remote directory
      f.cwd(f'/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.0.0/1km/{param}/mon/v20181126/')
      # define filename
      file=f'{param}_hadukgrid_uk_1km_mon_{year}01-{year}12.nc'
      # get the remote file to the local directory
      f.retrbinary("RETR %s" % file, open(file, "wb").write)

    # Close FTP connection
    f.close()
"""
function retrieve_HadUK(param::String, from_year::Int64, to_year::Int64, directory::String = ".")
    py"retrieve_HadUK"(param, from_year, to_year, directory)
end
