import numpy as np
import pyLCIO
import ROOT
from math import *
import argparse
import os
import random
from plothelper import *
from datetime import datetime

# setup plotter
plt = PlotHelper()

sensorAngles = np.arange(-np.pi-np.pi/8,np.pi+2*np.pi/8,np.pi/8)

def getYlocalAndGamma(x,y):
    # Get exact angle gamma of the hit position
    gammaP=np.arctan2(y,x)

    # Get two sensor angles that are closest to the exact angle
    diff = np.abs(sensorAngles-gammaP)
    index1 = np.argmin(diff)
    gamma1=sensorAngles[index1]
    diff[index1]=3*np.pi
    index2 = np.argmin(diff)
    gamma2=sensorAngles[index2]

    # Rotate x coordinate of the point by each option for gamma
    x1=x*np.cos(-gamma1)-y*np.sin(-gamma1)
    y1=y*np.cos(-gamma1)+x*np.sin(-gamma1)
    x2=x*np.cos(-gamma2)-y*np.sin(-gamma2)
    y2=y*np.cos(-gamma2)+x*np.sin(-gamma2)

    # Determine which x is closest to expected value
    xTrue=30.16475324197002

    diff1=abs(x1-xTrue)
    diff2=abs(x2-xTrue)
    
    # If both x1 and x2 are really close to the ex
    if diff1 < 0.5 and diff2 < 0.5:
        if y1>8.5 or y1<-4.5:
            index=index2
        else:
            index=index1
            
    elif diff1<diff2:
        index=index1
    else:
        index=index2

    if index==index1:
        yentry=y1
    else:
        yentry=y2
    
    ylocal=-round(yentry/25e-3)*25e-3
    # at some point, add limits to possible ROIs

    gamma=sensorAngles[index]
    # Shift range of gamma to 0 to 2 pi
    if gamma<0:
        gamma+=2*np.pi
    
    return ylocal, gamma

# user options
parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input_file", help="Input file", type=str)
parser.add_argument("-odir", "--output_folder", help="Output folder", type=str)
parser.add_argument("-f", "--float_precision", help="Floating point precision", default=5, type=int)
parser.add_argument("-t", "--track_total", help="Total number of tracks to simulate (for BIB and signal individually)", default=100, type=int)
parser.add_argument("-b", "--bin_size", help="Number of tracks per tracklist", default=1, type=int) 
parser.add_argument("-p", "--plot", help="Include if you want to make plots at this stage", action='store_true')
parser.add_argument("-flp", "--flp", help="Direction of sensor (1 for FE side out, 0 for FE side down)", default=0, type=int)
parser.add_argument("-sig", "--signal", help="Are you generating signal?", action='store_true')
parser.add_argument("-bmm", "--bib_mm", help="Are you running mm bib?", action='store_true')
parser.add_argument("-bmp", "--bib_mp", help="Are you running mp bib?", action='store_true')

ops = parser.parse_args()

# get a list of files based on what you want to make tracklists for
if ops.signal:
    if not ops.input_file: 
        raise Exception("You must specify an input file for signal")
    if not ops.input_file.endswith(".slcio"):
        raise Exception("Input file must be a .slcio file")
    file_list=[ops.input_file]
    output_file_form="signal_tracks_*.txt"
elif ops.bib_mm or ops.bib_mp:
    if ops.bib_mp:
        directory_path = "/cvmfs/public-uc.osgstorage.org/ospool/uc-shared/public/futurecolliders/BIB10TeV/sim_mp_pruned/" 
        output_file_form="bib_mp_tracks_*.txt"
    else:
        directory_path = "/cvmfs/public-uc.osgstorage.org/ospool/uc-shared/public/futurecolliders/BIB10TeV/sim_mm_pruned/" 
        output_file_form="bib_mm_tracks_*.txt"
    file_list=os.listdir(directory_path)
    file_list = [os.path.join(directory_path, file) for file in file_list]
else: 
    raise Exception("You must include one of the flags for signal (-sig), bib mm (-bmm), or bib mp (-bmp)")

output_file=f"{ops.output_folder}/{output_file_form}"
    
# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
ROOT.gROOT.SetBatch()

# check if output directory exists
#if not os.path.isdir(ops.output_file):
#    raise Exception(f"Directory {ops.output_file} does not exist")

# convert this to writing to a log file
#print(f"Getting tracks from file: {ops.input_file}\n")
#print(f"Allowed PIDs: {ops.allowedPIDS}\n")

tracks=[]

track_count=0
break_loop=False
count = 0
# BIB processing can skip into the large CVMFS file lists. Signal mode has a
# single detector-sim file, so starting at 440 silently drops all signal hits.
fileCountStart = 0 if ops.signal else 440
fileCountLimit = 1000
for file_path in file_list:
    if break_loop:
            break
    count +=1
    if count<fileCountStart:
        continue
    print("count", count, " ", datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    if count>fileCountLimit:
        print("count>fileCountLimit, so breaking")
        break
    # Get the full path to the file
    if not file_path.endswith(".slcio"):
        continue
    print("Processing file: ", file_path)
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file_path)

    for ievt,event in enumerate(reader):
        
        print("Processing event %i."%ievt)

        # Print all the collection names in the event
        collection_names = event.getCollectionNames()

        # Get vertex barrel hits
        vtxBarrelHits = event.getCollection("VertexBarrelCollection")

        for hit in vtxBarrelHits: 
            if track_count>ops.track_total-1:
                print(f"\nReached total number of tracks ({ops.track_total}). Breaking loop.\n")
                break_loop=True
                break

            position = hit.getPosition()
            
            x,y,hit_z = position[0],position[1],position[2]
            #rxy=(x**2+y**2)**0.5
            t=hit.getTime()

            # Get layer
            encoding = vtxBarrelHits.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
            decoder = pyLCIO.UTIL.BitField64(encoding)
            cellID = int(hit.getCellID0())
            decoder.setValue(cellID)
            #detector = decoder["system"].value()
            layer = decoder['layer'].value()
            #side = decoder["side"].value()

            if layer!=0: continue # first layer only

            if ops.signal:
                # get the particle that caused the hit
                mcp = hit.getMCParticle()
                hit_pdg = mcp.getPDG() if mcp else None
                hit_id = mcp.id() if mcp else None

                if hit_id is None:
                    print("No MCParticle associated with hit, skipping.")
                    continue

                if abs(hit_pdg) != 13: continue

                # momentum at production
                mcp_p = mcp.getMomentum()
                mcp_tlv = ROOT.TLorentzVector()
                mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())


                # momentum at hit
                hit_p = hit.getMomentum()
                hit_tlv = ROOT.TLorentzVector()
                hit_tlv.SetPxPyPzE( hit_p[0], hit_p[1], hit_p[2], mcp.getEnergy())
            else:
                #momentum at hit
                m_electron=.000511 #Gev
                hit_p = hit.getMomentum() #Gev
                hit_tlv = ROOT.TLorentzVector()
                hit_tlv.SetPxPyPzE( hit_p[0], hit_p[1], hit_p[2], np.sqrt(m_electron**2+hit_tlv.P()**2))

                hit_pdg=11 #random.choice(11,-11)?

                if True:
                    p1Calc = np.sqrt(hit_p[2]*hit_p[2]+hit_tlv.Pt()*hit_tlv.Pt())
                    if np.abs(p1Calc/hit_tlv.P() -1)>0.000001:
                        print(f"hit momentum {hit_p[2]}, pt {hit_tlv.Pt()}, p {hit_tlv.P()}, sqrt(pt^2 + pz^2 {p1Calc}")

            if ops.plot:
                # plt.plot1D("hit_mcp_e"  ,";mcp e [GeV];hits" , mcp_tlv.E(), 100, 0, 0.2)
                # plt.plot1D("hit_mcp_pt"  ,";mcp pt [GeV];hits" , mcp_tlv.Pt(), 100, 0, 0.2)
                # plt.plot1D("hit_mcp_eta" ,";mcp eta;hits" , mcp_tlv.Eta(), 100, -3.2, 3.2)
                # plt.plot1D("hit_mcp_theta" ,";mcp theta;hits" , mcp_tlv.Theta(), 100, 0, 3.2)
                # plt.plot1D("hit_mcp_phi" ,";mcp phi;hits" , mcp_tlv.Phi(), 100, -3.2, 3.2)
                plt.plot1D("hit_pz"  ,";incident pz [GeV];hits" , hit_p[2], 100, 0, 0.2)
                plt.plot1D("hit_pz_Long"  ,";incident pz [GeV];hits" , hit_p[2], 100, 0, 1)
                plt.plot1D("hit_p"  ,";incident p [GeV];hits" , hit_tlv.P(), 100, 0, 0.2)
                plt.plot1D("hit_p_Long"  ,";incident p [GeV];hits" , hit_tlv.P(), 100, 0, 1)
                plt.plot1D("hit_pt"  ,";incident pt [GeV];hits" , hit_tlv.Pt(), 100, 0, 0.2)
                plt.plot1D("hit_pt_Long"  ,";incident pt [GeV];hits" , hit_tlv.Pt(), 100, 0, 1)
                plt.plot1D("hit_eta" ,";incident eta;hits" , hit_tlv.Eta(), 100, -3.2, 3.2)
                plt.plot1D("hit_theta" ,";incident theta;hits" , hit_tlv.Theta(), 100, 0,3.2)
                plt.plot1D("hit_phi" ,";incident phi;hits" , hit_tlv.Phi(), 100, -3.2, 3.2)
                plt.plot1D("hit_eDep" ,";incident e deposit [MeV];hits" , hit.getEDep()*1000, 100, 0, 0.5)
                # plt.plot1D("hit_mcp_prodrxy" ,";mcp prod rxy [mm];hits" , prodrxy, 100, 0,150)
                # plt.plot1D("hit_mcp_prodz"   ,";mcp prod z [mm];hits" , prodz, 100, -1000,1000)
                # plt.plot1D("hit_mcp_endrxy"  ,";mcp end rxy [mm];hits" , endrxy, 100, 0, 150)
                # plt.plot1D("hit_mcp_endz"    ,";mcp end z [mm];hits" , endz, 100, -1000, 1000)

                plt.plot1D("hit_time"    ,";hit time [ns];hits" , t, 100, -1, 5)
                # print("phi particle,hit {:.2f} {:.2f}".format(hit_tlv.Phi(), mcp_tlv.Phi()))
                # print("eta particle,hit {:.2f} {:.2f}".format(hit_tlv.Eta(), mcp_tlv.Eta()))

                # double check if any bugs
                phi = hit_tlv.Phi()
                theta = hit_tlv.Theta()

                plt.plot1D("hit_phi"    ,";cota;hits" , phi, 100, -10,10)
                plt.plot1D("hit_theta"    ,";cotb;hits" , theta, 100, -10,10)
                plt.plot1D("hit_t"    ,";t;hits" , t, 100, -1,10)

            p = hit_tlv.P()
            pt = hit_tlv.Pt()
            
            ylocal, gamma0 = getYlocalAndGamma(x,y)
            zglobal = round(hit_z/25e-3)*25e-3 # round to nearest pixel
            
            # Alternative cota and cotb calculation
            cota = hit_p[2]/(hit_p[0]*cos(gamma0)+hit_p[1]*sin(gamma0))
            cotb = (hit_p[0]*sin(gamma0)-hit_p[1]*cos(gamma0))/(hit_p[0]*cos(gamma0)+hit_p[1]*sin(gamma0))

            # If we are unflipped, we must adjust alpha and beta, and flip y-local
            if ops.flp == 0:
                ylocal *= -1
                cota *= -1

            # Skip tracks with momentum 0
            if round(p, ops.float_precision)==0 or round(pt, ops.float_precision)==0:
                continue
            
            if np.abs(p) < 0.00005:  
                print("bad momentum!!")
                continue; #added because pixelav can't handle these.
            elif np.abs(pt) < 0.00005: 
                print("bad pt momentum!!")
                continue; #added because pixelav can't handle these.
            else:
            
                track = [cota, cotb, p, ops.flp, ylocal, zglobal, pt, t, hit_pdg]
                tracks.append(track)
                track_count+=1

binsize = ops.bin_size
float_precision = ops.float_precision
numFiles = int(np.ceil(len(tracks)/binsize))

print(f"Writing {len(tracks)} tracks to {numFiles} files with {binsize} tracks per file\n")

for fileNum in range(numFiles):
    with open(output_file.replace("*", str(fileNum)), 'w') as file:
        for track in tracks[fileNum*binsize:(fileNum+1)*binsize]:

            # set flp to an int
            track = list(track)

            formatted_sublist = [f"{element:.{float_precision}f}" if isinstance(element, float) else element for element in track]
            line = ' '.join(map(str, formatted_sublist)) + '\n'
            file.write(line)
