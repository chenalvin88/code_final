import numpy as np
from datahandling import filehandling
import matplotlib.pyplot as plt
import matplotlib.image as img
import cv2
from tqdm import tqdm as tqdm

data = filehandling("500_e7.csv")
a = data.extract('plateifu')
b = data.extract('plateifu',selection='hasjpa')
name_list = [e for e in tqdm(data.extract('plateifu')) if e not in data.extract('plateifu',selection='hasjpa')]
# ra_list = data.extract('objra')
# dec_list = data.extract('objdec')
file_out = open("download.sh","w")

# for (plateifu,ra,dec) in zip(name_list,ra_list,dec_list):
# for (plateifu,ra,dec) in zip(name_list[0],ra_list[0],dec_list[0]):
    # print(name_list[-4:])
for i,plateifu in enumerate(name_list[2:]):
    temp = plateifu.replace('-', ' ').split(' ')
    plate = int(temp[0])
    ifu = int(temp[1])
    # ra = data[(plateifu,'objra','float')]
    # dec = data[(plateifu,'objdec','float')]

    # manga dap
    # download_path1 = f"/Volumes/SDrive/yenting_pa_alignment/MANGA/SPX-MILESHC-MASTARSSP/manga-{plateifu}-MAPS-SPX-MILESHC-MASTARSSP.fits.gz"
    # download_path1 = f"/Volumes/SDrive/yenting_pa_alignment/MANGA/VOR10-MILESHC-MASTARSSP/manga-{plateifu}-MAPS-VOR10-MILESHC-MASTARSSP.fits.gz"
    download_path1 = f"/Volumes/SDrive/yenting_pa_alignment/MANGA/HYB10-MILESHC-MASTARSSP/manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz"
    # download_path1 = f"/Volumes/SDrive/yenting_pa_alignment/MANGA/HYB10-MILESHC-MASTARHC2/manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz"
    # download_url1 = f"--user='sdss' --password='2.5-meters' https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-11/SPX-MILESHC-MASTARSSP/{plate}/{ifu}/manga-{plateifu}-MAPS-SPX-MILESHC-MASTARSSP.fits.gz"
    # download_url1 = f"--user='sdss' --password='2.5-meters' https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-11/VOR10-MILESHC-MASTARSSP/{plate}/{ifu}/manga-{plateifu}-MAPS-VOR10-MILESHC-MASTARSSP.fits.gz"
    download_url1 = f"--user='sdss' --password='2.5-meters' https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-11/HYB10-MILESHC-MASTARSSP/{plate}/{ifu}/manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz"
    # download_url1 = f"--user='sdss' --password='2.5-meters' https://data.sdss.org/sas/mangawork/manga/spectro/analysis/MPL-11/HYB10-MILESHC-MASTARHC2/{plate}/{ifu}/manga-{plateifu}-MAPS-HYB10-MILESHC-MASTARHC2.fits.gz"
    file_out.write(f"wget -O {download_path1} {download_url1}\n")
    
    # manga pipe 3dhttps://data.sdss.org/sas/mangawork/manga/sandbox/pipe3d/v3_1_1/3.1.1/10001/manga-10001-12701.Pipe3D.cube.fits.gz
    # download_path1 = f"/Volumes/SDrive/yenting_pa_alignment/MANGA/Pipe3D/manga-{plateifu}.Pipe3D.cube.fits.gz"
    # download_url1 = f"--user='sdss' --password='2.5-meters' https://data.sdss.org/sas/mangawork/manga/sandbox/pipe3d/v3_1_1/3.1.1/{plate}/manga-{plateifu}.Pipe3D.cube.fits.gz"
    # file_out.write(f"wget -O {download_path1} {download_url1}\n")

    # # nvss
    # download_path2 = f"/Volumes/SDrive/yenting_pa_alignment/sauron_and_vla/nvss/name={plateifu}_ra={ra:.3f}_dec={dec:.3f}_nvss.jpg"
    # download_path2 = f"/Volumes/SDrive/yenting_pa_alignment/vla_proposal/selected/nvss/name={plateifu}_nvss.jpg"
    # download_url2 = f"http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?position={plateifu}&survey=nvss&pixels=20%2C20&sampler=Clip&scaling=Sqrt&size=default&projection=Tan&coordinates=J2000.0&return=fits"
    # download_url2 = f"http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?position={ra:.5f}%%2C{dec:.5f}&survey=nvss&pixels=60%%2C60&sampler=Clip&scaling=Sqrt&size=default&projection=Tan&coordinates=J2000.0&return=fits"
    # file_out.write(f"wget -O {download_path2} '{download_url2}'\n")
    
    # # first
    # download_path3 = f"/Volumes/SDrive/yenting_pa_alignment/sauron_and_vla/first/name={plateifu}_ra={ra:.3f}_dec={dec:.3f}_first.jpg"
    # download_path3 = f"/Volumes/SDrive/yenting_pa_alignment/vla_proposal/selected/first/name={plateifu}_first_b.jpg"
    # download_url3 = f"http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?position={plateifu}&survey=first&pixels=167%2C167&sampler=Clip&scaling=Sqrt&size=default&projection=Tan&coordinates=J2000.0&return=fits" 
    # download_url3 = f"https://third.ucllnl.org/cgi-bin/firstimage?RA={ra}&Dec={dec}&Equinox=J2000&ImageSize=4.5&MaxInt=10&GIF=1&Download=1" # unit of image size is arcmin, RA (hours) DEC (deg)
    # download_url3 = f"http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?position={ra:.5f}%%2C{dec:.5f}&survey=first&pixels=60%%2C60&sampler=Clip&scaling=Sqrt&size=default&projection=Tan&coordinates=J2000.0&return=fits"
    # file_out.write(f"wget -O {download_path3} '{download_url3}'\n")
    
    # # vlass
    # download_path4 = f"/Volumes/SDrive/yenting_pa_alignment/vla_proposal/selected/vlass/name={plateifu}_vlass.jpg"
    # download_url4 = f"http://legacysurvey.org/viewer/fits-cutout?ra={ra}&dec={dec}&size=300&layer=vlass"
    # file_out.write(f"wget -O {download_path4} '{download_url4}'\n")

    # # sdss
    # download_path5 = f"/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/sdss/name={plateifu}_sdss.jpg"
    # download_path5 = f"/Volumes/SDrive/yenting_pa_alignment/vla_proposal/massive/sdss/name={plateifu}_sdss.jpg"
    # download_url5 = f"http://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?position={plateifu}&survey=sdssr&pixscale=0.396&pixels=758%2C758&sampler=Clip&scaling=Sqrt&size=default&projection=Tan&coordinates=J2000.0&return=jpeg"
    # download_url5 = f"http://legacysurvey.org/viewer/fits-cutout?ra={ra:.5f}&dec={dec:.5f}&pixscale=0.396&size=60&layer=sdss&bands=r"
    # file_out.write(f"wget -O {download_path5} '{download_url5}'\n")

    # # optical image
    # download_path6 = f"/Volumes/SDrive/yenting_pa_alignment/vla_proposal/selected/optical/name={plateifu}_optical.jpg"
    # download_path6 = f"/Volumes/SDrive/yenting_pa_alignment/radio_and_optical_images/optical2/name={plateifu}_optical.jpg"
    # download_url6 = f"http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx\?ra={ra:.5f}&dec={dec:.5f}&scale=0.396&width=758&height=758"
    # opt explanation : http://skyserver.sdss.org/dr1/en/help/docs/api.asp
    # download_url6 = f"http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx\?ra={ra:.5f}&dec={dec:.5f}&scale=0.396&width=273&height=273&opt=G"
    # http://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?ra=184.9511&dec=-0.8754&scale=0.4&height=512&width=512&opt=GO 
    # file_out.write(f"wget -O {download_path6} '{download_url6}'\n")
                    
    # L=9
    # if i==0:
        # plt.figure(figsize=(7,79/7*9/8),constrained_layout=False)
        # plt.subplots_adjust(hspace=0.0,wspace=0.0,left=0.2,right=0.8)
    # plt.subplot(1,4,2).imshow(img.imread(download_path2))
    # plt.subplot(L,4,2+i*4).imshow(np.array(cv2.cvtColor(cv2.imread(download_path2),cv2.COLOR_BGR2RGB)))
    # plt.title('nvss')
    # plt.axis('off')

    # plt.subplot(1,4,3).imshow(img.imread(download_path3))
    # plt.subplot(L,4,3+i*4).imshow(np.array(cv2.cvtColor(cv2.imread(download_path3),cv2.COLOR_BGR2RGB)))
    # plt.title('first')
    # plt.axis('off')
    # try:
        # plt.subplot(L,4,4+i*4).imshow(img.imread(download_path4))
        # plt.title('vlass1.2')
    # except  Exception:
        # pass
    # plt.axis('off')
    # try:
    #     plt.subplot(234).imshow(img.imread(download_path5))
    #     plt.title('sdss_r')
    # except Exception:
    #     pass
    # plt.axis('off')
    # plt.subplot(L,4,1+i*4).imshow(img.imread(download_path6))
    # plt.title('sdss')
    # plt.annotate(plateifu, (0,100),color='white')
    # plt.axis('off')
    # plt.tight_layout()
    # if i==L-1:
        # plt.savefig(f'/Volumes/SDrive/yenting_pa_alignment/vla_proposal/scientific_justification/tight2.png',bbox_inches='tight')
        # plt.show()
    # plt.close()


file_out.close()