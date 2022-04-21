# Importing modules
import sys
import os
sys.path.append(os.path.abspath(R"C:\Users\SoPa\Documents\GitHub\spectrum_awg\RPB"))
from pyspcm import *
from spcm_tools import *
import ctypes
import numpy as np
import numexpr as ne
from multiprocessing import Process, Value, Array
from timeit import default_timer as timer

def get_freq_and_amp(sign,frequencies,amplitudes):
    """
    This function calculates the frequencies for the immediately next step.
    """
    freq_update_step = 50e3 # in Hz. frequency step to take while moving the tweezers.
    freq = frequencies+(sign*freq_update_step)
    amp = amplitudes
    return freq, amp

def calculate_target_frequency(start_frequencies,tweezer_occu):
    """
    This function calculates the target frequency arrangement. Target is
    same for all intil configs with same number of tweezers and loaded
    atoms. The target arrangement is assembled around the center of full array.
    """
    inds,=np.asarray(tweezer_occu==1).nonzero() #
    center=np.take(inds, inds.size//2)
    z,=np.asarray(inds==center).nonzero() #
    new_inds=np.arange(center-z[0],center-z[0]+inds.size)
    new_inds=new_inds+(tweezer_occu.size//2)-center
    return start_frequencies[new_inds]

def calculate_signal(t, frequencies,amplitudes):
    """
    This function calculates the total data using numpy
    vectorization and broadcasting.
    """
    return np.sum(amplitudes[:,None]*np.sin(2*np.pi*frequencies[:,None]*t[None,:]), axis=0).astype('int16')

def calculate_signal_fast(t,frequencies,amplitudes):
    """
    This function calculates the total data using numpy
    and numexpr. Performance is better than numpy alone.
    """
    return np.sum(amplitudes[:,None]*ne.evaluate('sin(arg)',{'arg':2*np.pi*frequencies[:,None]*t[None,:]}), axis=0).astype('int16')

def jump_between(check, free, data_1_arr, data_2_arr, t):
    """
    This function runs in an independent side process.
    It manipulates the shared arrays and booleans to stream properly to AWG.
    """
    start_frequencies = np.arange(42., 62.0, 1.0)*1e6 # frequencies for all tweezers produced intitally
    start_amplitudes = np.linspace(10000//start_frequencies.size, 10000//start_frequencies.size, start_frequencies.size)
    target_frequencies = np.ones(start_frequencies.size) # frequencies of the defect free array for a particular load configuration
    frequencies = np.ones(start_frequencies.size) # actual frequencies that are to be streamed to the AWG
    amplitudes = np.ones(start_frequencies.size)
    sign = np.ones(start_frequencies.size)
    file_read=False # local boolean variable for tracking the occupancy text file
    no_atoms=False # local boolean variable for representing an empty array
    time_once=True # local boolean variable for measuring sorting performance
    time_diff=0.0 # local float variable for measuring sorting performance
    textfile_path=R'C:\Users\SoPa\Documents\tweezer_occupancy.txt'
    while True:
        if file_read==False:
            try:
                with open(textfile_path, "r") as text_file:
                    # time_diff = timer()
                    tweezer_occu=text_file.read()
                file_read=True
                tweezer_occu=np.asarray(list(tweezer_occu), dtype=int) # read the tweezer occupancy text file as a numpy array
                assert tweezer_occu.size == start_frequencies.size, 'Tweezer occupancy file invalid'
                if np.any(tweezer_occu):
                    frequencies,amplitudes = start_frequencies[tweezer_occu==1],start_amplitudes[tweezer_occu==1]
                    target_frequencies = calculate_target_frequency(start_frequencies,tweezer_occu)
                else:
                    no_atoms=True
            except Exception as e:
                print(e)
                pass
            continue
        if no_atoms:
            continue
        if free.value==1:
            free.value=0
            freq_diff=target_frequencies-frequencies
            sign=np.sign(freq_diff) # array to store whether each frequency os above or below the target
            if np.any(sign): # this executes until target frequencies are reached
                frequencies,amplitudes = get_freq_and_amp(sign,frequencies,amplitudes)
            else: # this executes when target frequencies are reached
                if time_once: # this executes once, when target frequencies are reached
                    # time_diff=timer()-time_diff
                    # print(time_diff)
                    time_once=False
                    pass
            if check.value == 0:
                try:
                    temp_arr_1=calculate_signal_fast(t, frequencies,amplitudes)
                    ctypes.memmove(ctypes.byref(data_1_arr.get_obj()), temp_arr_1.ctypes.data, temp_arr_1.nbytes)
                except Exception as e:
                    print(e)
                check.value = 1
            else:
                try:
                    temp_arr_2=calculate_signal_fast(t, frequencies,amplitudes)
                    ctypes.memmove(ctypes.byref(data_2_arr.get_obj()), temp_arr_2.ctypes.data, temp_arr_2.nbytes)
                except Exception as e:
                    print(e)
                check.value = 0

# Buffer data transfer function
def vCalcNewData (pnBuffer, lNumCh, llSamplePos, llNumSamples,final_data):
    lStartPosInBuffer_bytes = (llSamplePos % lPreCalcLen) * 2 * lNumCh #why is the factor 2 needed here? int16 2 byte long
    lToCopy_bytes = llNumSamples * 2 * lNumCh
    lPreCalcLen_bytes = lPreCalcLen * 2 * lNumCh
    lAlreadyCopied_bytes = 0
    while lAlreadyCopied_bytes < lToCopy_bytes:
        # copy at most the pre-calculated data
        lCopy_bytes = lToCopy_bytes - lAlreadyCopied_bytes
        if lCopy_bytes > lPreCalcLen_bytes - lStartPosInBuffer_bytes:
            lCopy_bytes = lPreCalcLen_bytes - lStartPosInBuffer_bytes

        # copy data from pre-calculated buffer to DMA buffer. The get_obj() function is for shared Array variable.
        ctypes.memmove (cast (pnBuffer, c_void_p).value + lAlreadyCopied_bytes, cast (final_data.get_obj(), c_void_p).value + lStartPosInBuffer_bytes, lCopy_bytes)
        lAlreadyCopied_bytes += lCopy_bytes
        lStartPosInBuffer_bytes = 0

if __name__ == '__main__':
    # Initiate the card
    spcm_vClose (-1) #this line closes the card first for the case that an error occured which stopped the compiling before the card is closed.
    # open card
    # uncomment the second line and replace the IP address to use remote cards like in a generatorNETBOX
    hCard = spcm_hOpen (create_string_buffer (b'/dev/spcm0'))
    #hCard = spcm_hOpen (create_string_buffer (b'TCPIP::192.168.1.10::inst0::INSTR'))
    if hCard == None:
        sys.stdout.write("no card found...\n")
        exit ()
    # read type, function and SN and check for D/A card
    lCardType = int32 (0)
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP, byref (lCardType))
    lSerialNumber = int32 (0)
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO, byref (lSerialNumber))
    lFncType = int32 (0)
    spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE, byref (lFncType))
    sCardName = szTypeToName (lCardType.value)
    if lFncType.value == SPCM_TYPE_AO:
        sys.stdout.write("Found: {0} sn {1:05d}\n".format(sCardName,lSerialNumber.value))
    else:
        sys.stdout.write("This is an example for D/A cards.\nCard: {0} sn {1:05d} not supported by example\n".format(sCardName,lSerialNumber.value))
        spcm_vClose (hCard);
        exit ()

    # Setup the card
    lPreCalcLen = int(0) # in samples
    Rate = 500*(10**6) #in Hz, 1250MHz is max for this card
    # set samplerate to Rate, no clock output
    spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, Rate)
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,   0)
    # driver might have adjusted the sampling rate to the best-matching value, so we work with that value
    SRate = int64 (0)
    spcm_dwGetParam_i64 (hCard, SPC_SAMPLERATE, byref (SRate))

    # setup the mode
    qwChEnable = uint64 (CHANNEL0)
    spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,    SPC_REP_FIFO_SINGLE)#SPC_REP_FIFO_SINGLE)
    spcm_dwSetParam_i64 (hCard, SPC_CHENABLE,    qwChEnable)
    spcm_dwSetParam_i64 (hCard, SPC_SEGMENTSIZE, 4096) # used to limit amount of replayed data if SPC_LOOPS != 0
    spcm_dwSetParam_i64 (hCard, SPC_LOOPS,       0) # continuous replay
    lSetChannels = int32 (0)
    spcm_dwGetParam_i32 (hCard, SPC_CHCOUNT, byref (lSetChannels))
    lBytesPerSample = int32 (0)
    spcm_dwGetParam_i32 (hCard, SPC_MIINST_BYTESPERSAMPLE,  byref (lBytesPerSample))

    # setup the trigger mode (SW trigger, no output)
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,      SPC_TMASK_SOFTWARE)
    #spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,      SPC_TMASK_NONE)
    #spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK, SPC_TMASK_EXT0)
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK, 0)
    #spcm_dwSetParam_i32 (hCard, SPC_TRIG_EXT0_LEVEL0, 1000)
    #spcm_dwSetParam_i32 (hCard, SPC_TRIG_EXT0_MODE, SPC_TM_HIGH)
    #spcm_dwSetParam_i32 (hCard, SPC_CH0_STOPLEVEL, SPCM_STOPLVL_HIGH)
    ########
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ORMASK0,  0)
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ORMASK1,  0)
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ANDMASK0, 0)
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ANDMASK1, 0)
    spcm_dwSetParam_i32 (hCard, SPC_TRIGGEROUT,       0)

    # setup all channels
    for i in range (0, lSetChannels.value):
        spcm_dwSetParam_i32 (hCard, SPC_AMP0 + i * (SPC_AMP1 - SPC_AMP0), int32 (1000))
        spcm_dwSetParam_i32 (hCard, SPC_ENABLEOUT0 + i * (SPC_ENABLEOUT1 - SPC_ENABLEOUT0),  int32(1))

    # Data calculation
    ne.set_num_threads(8) #number of threads for numexpr
    frequencies =  np.arange(42., 62.0, 1.0)*1e6 # frequencies for all tweezers produced intitally
    amplitudes = np.linspace(10000//frequencies.size, 10000//frequencies.size, frequencies.size)
    freq_resolution=50e3 # frequency resolution desired
    end_time = 1/freq_resolution # time required to achieve desired frequency resolution
    number_of_samples = int(SRate.value*end_time)
    t=np.linspace(0.0, end_time, num=number_of_samples, endpoint=False) # time points array
    data_1=calculate_signal_fast(t, frequencies,amplitudes) # data points array
    data_2=calculate_signal_fast(t, frequencies,amplitudes) # data points array
    lPreCalcLen = len(data_1)
    print('data size in MB', lPreCalcLen*2/1024/1024)

    # setup hardware buffer (card memory)
    llHWBufSize= uint64 (32*1024*1024) # Do not make too big to reduce latency
    spcm_dwSetParam_i64 (hCard, SPC_DATA_OUTBUFSIZE, llHWBufSize);
    spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_WRITESETUP);

    # setup software buffer or PC memory
    lNotifySize_bytes = int32(512*1024) # 1 MByte
    qwBufferSize = uint64 (30*1024*1024) # For simplicity qwBufferSize should be a multiple of lNotifySize_bytes

    # we try to use continuous memory if available and big enough
    pvBuffer = c_void_p ()
    qwContBufLen = uint64 (0)
    spcm_dwGetContBuf_i64 (hCard, SPCM_BUF_DATA, byref(pvBuffer), byref(qwContBufLen))
    sys.stdout.write ("ContBuf length: {0:d}\n".format(qwContBufLen.value))
    if qwContBufLen.value >= qwBufferSize.value:
        sys.stdout.write("Using continuous buffer\n")
    else:
        pvBuffer = pvAllocMemPageAligned (qwBufferSize.value)
        sys.stdout.write("Using buffer allocated by user program\n")

    # we calculate data for all enabled channels, starting at sample position 0, and fill the complete DMA buffer
    check = Value('b', True) # synchronized shared Boolean variable
    free = Value('b', False) # synchronized shared Boolean variable
    data_1_arr = Array('i', range(lPreCalcLen)) # synchronized shared Array variable
    data_2_arr = Array('i', range(lPreCalcLen)) # synchronized shared Array variable
    ctypes.memmove(ctypes.byref(data_1_arr.get_obj()), data_1.ctypes.data, data_1.nbytes)
    ctypes.memmove(ctypes.byref(data_2_arr.get_obj()), data_2.ctypes.data, data_2.nbytes)
    qwSamplePos = 0
    lNumAvailSamples = (qwBufferSize.value // lSetChannels.value) // lBytesPerSample.value
    vCalcNewData (pvBuffer, lSetChannels.value, qwSamplePos, lNumAvailSamples, data_1_arr)
    qwSamplePos += lNumAvailSamples

    # we define the buffer for transfer and start the DMA transfer
    sys.stdout.write("Starting the DMA transfer and waiting until data is in board memory\n")
    spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, lNotifySize_bytes, pvBuffer, uint64 (0), qwBufferSize)
    spcm_dwSetParam_i32 (hCard, SPC_DATA_AVAIL_CARD_LEN, qwBufferSize)
    spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA)

    # we'll start and wait until the card has finished or until a timeout occurs
    lStatus = int32(0)
    lAvailUser_bytes = int32(0)
    lPCPos = int32(0)
    lFillsize = int32(0)
    bStarted = False
    # start the side proces to run the function jump_between()
    p = Process(target=jump_between, args=(check, free, data_1_arr, data_2_arr, t))
    p.start()
    # continuously stream to AWG
    while True:
        dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA)
        if dwError != ERR_OK:
            if dwError == ERR_TIMEOUT:
                sys.stdout.write ("... Timeout\n")
            else:
                sys.stdout.write ("... Error: {0:d}\n".format(dwError))
                break;
        else:
            # start the card if the onboard buffer has been filled completely
            spcm_dwGetParam_i32 (hCard, SPC_FILLSIZEPROMILLE, byref (lFillsize));
            if lFillsize.value == 1000 and bStarted == False:
                sys.stdout.write("... data has been transferred to board memory\n")
                sys.stdout.write("\nStarting the card...\n")
                dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER)
                if dwError == ERR_TIMEOUT:
                    spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_STOP)
                    sys.stdout.write ("... Timeout at start\n")
                    break;
                bStarted = True
            else:
                sys.stdout.write ("... Fillsize: {0:d}/1000\n".format(lFillsize.value))
                pass
            spcm_dwGetParam_i32 (hCard, SPC_M2STATUS,            byref (lStatus))
            spcm_dwGetParam_i32 (hCard, SPC_DATA_AVAIL_USER_LEN, byref (lAvailUser_bytes))
            spcm_dwGetParam_i32 (hCard, SPC_DATA_AVAIL_USER_POS, byref (lPCPos))
            # calculate new data
            if lAvailUser_bytes.value >= lNotifySize_bytes.value:
                if check.value == 1:
                    free.value=1
                    # time_diff = timer()
                    pnData = (c_char * (qwBufferSize.value - lPCPos.value)).from_buffer (pvBuffer, lPCPos.value)
                    lNumAvailSamples = (lNotifySize_bytes.value // lSetChannels.value) // lBytesPerSample.value # to avoid problems with buffer wrap-arounds we fill only one notify size
                    vCalcNewData (pnData, lSetChannels.value, qwSamplePos, lNumAvailSamples,data_1_arr)
                else:
                    free.value=1
                    # time_diff = timer()
                    pnData = (c_char * (qwBufferSize.value - lPCPos.value)).from_buffer (pvBuffer, lPCPos.value)
                    lNumAvailSamples = (lNotifySize_bytes.value // lSetChannels.value) // lBytesPerSample.value # to avoid problems with buffer wrap-arounds we fill only one notify size
                    vCalcNewData (pnData, lSetChannels.value, qwSamplePos, lNumAvailSamples,data_2_arr)
                spcm_dwSetParam_i32 (hCard, SPC_DATA_AVAIL_CARD_LEN, lNotifySize_bytes)
                qwSamplePos += lNumAvailSamples
                # time_diff=timer()-time_diff
                # print(time_diff)
                free.value=0
    # Stop the card
    # send the stop command
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA)
    spcm_vClose (hCard);
    # End the side process
    p.terminate()
    p.join()
