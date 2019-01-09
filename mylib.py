# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 00:42:40 2018

@author: Mattia
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def edfToDataFrame(f):
    # number of rows
    n = f.getNSamples()[0]
    # number of columns
    m = f.signals_in_file
    # init a matrix
    sigbufs = np.zeros([n, m])
    # fill the data
    for i in np.arange(m):
        sigbufs[:,i] = f.readSignal(i)
    # create pandas dataframe
    data = pd.DataFrame(sigbufs)
    # rename columns
    data.columns= f.getSignalLabels()

    return data


def print_some_stuff(f):

    print("edfsignals: %i" % f.signals_in_file)
    print("file duration: %i seconds" % f.file_duration)
    print("startdate: %i-%i-%i" % (f.getStartdatetime().day,f.getStartdatetime().month,f.getStartdatetime().year))
    print("starttime: %i:%02i:%02i" % (f.getStartdatetime().hour,f.getStartdatetime().minute,f.getStartdatetime().second))
    print("recording: %s" % f.getPatientAdditional())
    print("patientcode: %s" % f.getPatientCode())
    print("gender: %s" % f.getGender())
    print("birthdate: %s" % f.getBirthdate())
    print("patient_name: %s" % f.getPatientName())
    print("patient_additional: %s" % f.getPatientAdditional())
    print("admincode: %s" % f.getAdmincode())
    print("technician: %s" % f.getTechnician())
    print("equipment: %s" % f.getEquipment())
    print("recording_additional: %s" % f.getRecordingAdditional())
    print("datarecord duration: %f seconds" % f.getFileDuration())
    print("number of datarecords in the file: %i" % f.datarecords_in_file)
    print("number of annotations in the file: %i" % f.annotations_in_file)
    # choose a channel
    channel = 3
    print("\nsignal parameters for the %d.channel:\n\n" % channel)
    print("label: %s" % f.getLabel(channel))
    print("samples in file: %i" % f.getNSamples()[channel])
    print("physical maximum: %f" % f.getPhysicalMaximum(channel))
    print("physical minimum: %f" % f.getPhysicalMinimum(channel))
    print("digital maximum: %i" % f.getDigitalMaximum(channel))
    print("digital minimum: %i" % f.getDigitalMinimum(channel))
    print("physical dimension: %s" % f.getPhysicalDimension(channel))
    print("prefilter: %s" % f.getPrefilter(channel))
    print("transducer: %s" % f.getTransducer(channel))
    print("samplefrequency: %f" % f.getSampleFrequency(channel))

    annotations = f.readAnnotations()
    for n in np.arange(f.annotations_in_file):
        print("annotation: onset is %f    duration is %s    description is %s" % (annotations[0][n],annotations[1][n],annotations[2][n]))

    buf = f.readSignal(channel)
    n = 200
    print("\nread %i samples\n" % n)
    result = ""
    for i in np.arange(n):
        result += ("%.1f, " % buf[i])
    print(result)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def plot_analysis(densities, clustering_coeffs, average_path_lengths):
    f, axarr = plt.subplots(2, sharex=True, figsize = (9,9))
    axarr[0].plot(densities, clustering_coeffs[::-1], color = "green")
    axarr[0].set_title('Global Indices vs Graph Density')
    axarr[0].set_xlabel("density")
    axarr[0].set_ylabel("clustering coefficient")
    axarr[0].grid(linestyle = ":")
    axarr[1].plot(densities, average_path_lengths[::-1], color = "red")
    axarr[1].set_xlabel("density")
    axarr[1].set_ylabel("average path length")
    axarr[1].grid(linestyle = ":")
    plt.show()
    
def write_inputFileForMotifAnalysis(G, file):
    with open(file, 'a') as file:
        for source in range(G.vcount()):
            targets = G.get_adjlist()[source]
            for target in targets:
                if target != source:
                    file.write(str(source) + "  " + str(target) +"  " + str(1)+"\n")
    file.close()
    return
    