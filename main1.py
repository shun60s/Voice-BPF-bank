#coding:utf-8

#
# This shows BPF bank output envlop between 2000Hz and 4500Hz step by 5Hz and gray scale image file output
# 2000から4500Hzの間で5Hz毎のBPF出力波形の包絡線（絶対値）を画像データ（グレースケール）として作成する
#
# It can re-load of the pre-calculated BPF bank data before.
#

import sys
import os
import argparse
import numpy as np
from scipy import signal
from scipy.io.wavfile import read as wavread
from scipy.io.wavfile import write as wavwrite
from matplotlib import pyplot as plt

from iir1 import *
from BPF import *


# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  matplotlib  2.1.1
#  scipy 1.0.0

class Class_sub2(object):
    def __init__(self, path0,  frame_num=None, NFRAME=640, NSHIFT=320, \
    Start_Freq=2000, End_Freq=4500, Delta_Freq=5, Q0=40.0,\
    path_bpfwav='BPF_out/BPF_bank_out', result_dir='result', save_png=False, \
    spectrum_leak=True):
        # initalize
        try:
            sr, y = wavread(path0)
        except:
            print ('error: wavread ')
            sys.exit()
        else:
            self.yg= y / (2 ** 15)
            self.sr= sr
            print ('sampling rate ', sr)
        
        # instance HPF to reduce effect of pitch and F1
        hpf=Class_IIR1(fc=800.0, btype='high', n_order=4, sampling_rate=self.sr)  # -24dB/oct
        self.hpf_wav= hpf.filtering( self.yg)
        if 1:  # set 1, if hpf wav write out to use Vocal-Tube-Estimation
            f, _= os.path.splitext( path0) # os.path.basename(path0))
            path_hpf = f + '_hpf' + '.wav'
            wavwrite( path_hpf, self.sr , ( self.hpf_wav * 2 ** 15).astype(np.int16))
            print ('save hpf_wav into ', path_hpf)
        
        
        # BPF bank
        self.Start_Freq= Start_Freq
        self.End_Freq= End_Freq
        self.Delta_Freq= Delta_Freq
        self.bpf_fc_list= np.array( np.append( np.arange(Start_Freq, End_Freq, Delta_Freq), End_Freq) )
        self.num_bpf_bank= len( self.bpf_fc_list)
        print ('number of BPF', self.num_bpf_bank )
        self.Q0= Q0
        f, _= os.path.splitext( os.path.basename(path0))
        self.path_bpfwav = path_bpfwav + '_' + f + '.npz'
        
        # check if directory exists 
        dir0=os.path.dirname(self.path_bpfwav)
        if len(dir0) > 0 and (not os.path.isdir(dir0)):  # if use directory 
            os.mkdir(dir0)
        
        if not os.path.isfile(self.path_bpfwav): # compute newly
            self.BPF_bank_wav= np.zeros( (self.num_bpf_bank, len(self.hpf_wav)) )
            self.BPF_bank_wav_envlop= np.zeros( (self.num_bpf_bank, len(self.hpf_wav)) )
            self.BPF_bank()
            np.savez(self.path_bpfwav, bwav=self.BPF_bank_wav, bwavenv=self.BPF_bank_wav_envlop)
            print ('saved data to...', self.path_bpfwav)
        else: # load precomuted data
            xdata= np.load(self.path_bpfwav)
            self.BPF_bank_wav= xdata['bwav']
            self.BPF_bank_wav_envlop= xdata['bwavenv']
            print ('load computed data from ', self.path_bpfwav)
            print ('remove (or rename) it if you want to make precomputed data newly !')
        
        
        # show specified frame
        self.NFRAME=NFRAME     # 640 sr=16Khz 40mS  # 400 sr=16Khz 25mS 
        self.NSHIFT=NSHIFT     # 320 sr=16Khz 20mS  # 160 sr=16khz 10mS
        
        if frame_num is not None and frame_num >= 0:
            # compute start point and end point of current l-th frame
            self.frame_num= frame_num
            sp= self.NSHIFT * self.frame_num
            ep= sp + self.NFRAME
            if ep > len(self.yg):
                ep= len(self.yg)
            # result display on screen
            self.base_pathf= get_path_name2( result_dir, path0, frame_num)
            self.make_figure(sp,ep, save_png=save_png)
            
            # experiment 1
            if spectrum_leak:
                self.gen_sin(sp,ep)
            
            
        else: # make all figures at once and save in result directory
            count= int(((len(self.yg) - ( self.NFRAME - self.NSHIFT)) / self.NSHIFT))
            for i in range (count):
                self.frame_num= i
                print ('frame_num ', self.frame_num)
                sp= self.NSHIFT * self.frame_num
                ep= sp + self.NFRAME
                if ep > len(self.yg):
                    ep= len(self.yg)
                # save figure as a file
                pathf= get_path_name( result_dir, path0, self.frame_num)
                # make figure
                self.make_figure(sp,ep, pathf)
        

    def BPF_bank(self,):
        # BPF process
        for i in range(self.num_bpf_bank):
            #print ('loop no. ',i)
            sys.stdout.write("\r%d" % i)
            sys.stdout.flush()
            
            bpfx=Class_BPF(fc=self.bpf_fc_list[i], Q=self.Q0,\
             sampling_rate=self.sr)
            bpf_wav= bpfx.iir2(self.hpf_wav)
            self.BPF_bank_wav[i]= bpf_wav
            self.BPF_bank_wav_envlop[i]= np.abs(signal.hilbert(bpf_wav)) # get envlop
            
    """
     experiment : bpf again by pure sin waveform with same envelop
                   to see leak spectrum caused by the envelop
    """
    def gen_sin(self,sp,ep):
        gsin0 =np.sin(np.arange(self.BPF_bank_wav_envlop.shape[1]) * \
        (2.0 * np.pi) *  (self.freq_pk312 / self.sr))
        
        gsin1= gsin0 * self.BPF_bank_wav_envlop[self.pk312]
        print ( gsin1.shape)
        plt.figure()
        plt.subplot(211)
        plt.title('sin ' + str( int(self.freq_pk312)) + 'Hz waveform with envlop temporal change')
        plt.plot( gsin1)
        plt.subplot(212)
        plt.title('BPF output envlop temporal change')
        plt.plot( self.BPF_bank_wav_envlop[self.pk312])
        plt.tight_layout()
        plt.show()
        
        #
        f, ext = os.path.splitext(self.path_bpfwav)
        path_bpfwav_sin = f + '_' + str(self.pk312) + ext
        print (path_bpfwav_sin)
        
        if not os.path.isfile(path_bpfwav_sin): 
            # BPF process
            BPF_bank_wav_sin= np.zeros( (self.num_bpf_bank, self.BPF_bank_wav_envlop.shape[1]) )
            BPF_bank_wav_envlop_sin= np.zeros( (self.num_bpf_bank, self.BPF_bank_wav_envlop.shape[1]) )
            for i in range(self.num_bpf_bank):
                #print ('loop no. ',i)
                sys.stdout.write("\r%d" % i)
                sys.stdout.flush()
                
                bpfy=Class_BPF(fc=self.bpf_fc_list[i], Q=self.Q0, sampling_rate=self.sr)
                bpf_wav= bpfy.iir2(gsin1)
                BPF_bank_wav_sin[i]= bpf_wav
                BPF_bank_wav_envlop_sin[i]= np.abs(signal.hilbert(bpf_wav))  # get envlop
            
            np.savez(path_bpfwav_sin, bwav=BPF_bank_wav_sin, bwavenv=BPF_bank_wav_envlop_sin)
            print ('saved data to...', self.path_bpfwav)
        else:
            ydata= np.load(path_bpfwav_sin)
            BPF_bank_wav_sin= ydata['bwav']
            BPF_bank_wav_envlop_sin= ydata['bwavenv']
            print ('load computed data from ', path_bpfwav_sin)
            print ('remove (or rename) it if you want to make precomputed data newly !')
        
        delta_x= ep - sp
        fig_sin= np.zeros( (self.num_bpf_bank, delta_x))
        wav_raw_sin= np.zeros( (self.num_bpf_bank, delta_x))
        for i in range( self.num_bpf_bank):
            fig_sin[i]= BPF_bank_wav_envlop_sin[i][sp:ep].copy()
            wav_raw_sin[i]= BPF_bank_wav_sin[i][sp:ep].copy()
        
        f=  fig_sin / np.amax( self.fig) # amax(self.fig) is the reference value
        fig_unit_sin = np.uint8(np.around(   np.where(f > 1.0, 1.0, f) * 255))
        
        spt= self.NSHIFT * self.frame_num / self.sr * 1000
        ept= (self.NSHIFT * self.frame_num  + delta_x) / self.sr * 1000
        
        # show comparison original and pure sin waveform
        fig, (ax1,ax2) = plt.subplots(ncols=2)
        
        ax1.set_title( 'frame no. ' + str(self.frame_num) +' : ' + str(spt) + ' to ' + str(ept) + ' msec'  )
        ax1.imshow( self.fig_rgb, origin='lower', extent=[0, delta_x, self.Start_Freq, self.End_Freq ],\
        aspect=0.5) #1.0)
        
        #ax2.set_title( 'frame no. ' + str(self.frame_num) +' : ' + str(spt) + ' to ' + str(ept) + ' msec'  )
        ax2.set_title('   spectrum leak around' + str( int(self.freq_pk312)) +'Hz')
        ax2.imshow( fig_unit_sin, cmap='gray', origin='lower', \
        extent=[0, delta_x, self.Start_Freq, self.End_Freq ], aspect=0.5) #1.0)
        #plt.tight_layout()
        plt.show()
    
    
    def gray2rgb(self, in_fig, wav_raw, red_mark_threshold=250, show_number=1, \
    freq_rang=None, save_png=False):
        #
        # show BPF output envlop of which value is highest
        #
        # The points of which value is more than red_mark_threshold becomes red
        
        if show_number > 0:
            #
            if freq_rang is not None:
                freq_rang= np.clip(np.array( freq_rang), self.Start_Freq, self.End_Freq)
                f_bottom= int((freq_rang[0] - self.Start_Freq) / self.Delta_Freq)
                f_upper  = int((freq_rang[1] - self.Start_Freq) / self.Delta_Freq)
                index_max= np.argmax(in_fig[f_bottom:f_upper,:] )
                p312=[ (index_max % in_fig.shape[1]) + f_bottom] 
                self.pk312= int (index_max / in_fig.shape[1]) + f_bottom
            else:
                index_max= np.argmax(in_fig)
                p312=[ index_max % in_fig.shape[1] ]
                self.pk312= int (index_max / in_fig.shape[1])
            
            self.freq_pk312= self.pk312 * self.Delta_Freq +  self.Start_Freq
            
            plt.figure()
            plt.subplot(211)
            plt.title('frame no. ' + str(self.frame_num) + ' maximum BPF ' + str(int(self.freq_pk312)) + ' : envelop and waveform')
            plt.xlabel('time step')
            d0= np.abs(in_fig[self.pk312,:])
            f0= np.arange(len(d0))
            index_peak= signal.find_peaks_cwt( d0 , [ 2 ,3 ])
            plt.plot( f0, d0)
            plt.plot( f0[index_peak], d0[ index_peak ] , "ro")
            plt.grid()
            # draw raw wave
            plt.subplot(212)
            plt.plot( wav_raw[self.pk312,:] )
            #plt.plot( wav_raw[self.pk312,0:250] )
            plt.grid()
            
            plt.tight_layout()
            
            pathf= self.base_pathf + '_time_step_' + '.png'
            if save_png:
                plt.savefig(pathf)
                print ( pathf)
            else:
                plt.show()
        
        # show BPF output
        for i in range (show_number):
            d0= np.abs(in_fig[ :, p312[0] + i])
            f0= np.arange(len(d0)) * self.Delta_Freq  + self.Start_Freq
            plt.figure()
            plt.title('frame no. ' + str(self.frame_num) + ' maximum BPF bank output:  '  + str(p312[0]) + '-' + str(i) )
            plt.xlabel('Hz')
            index_peak= signal.find_peaks_cwt( d0 , [ 2 ,3 ])
            plt.plot( f0, d0 )
            plt.plot( f0[index_peak], d0[ index_peak ] , "ro")
            
            pathf= self.base_pathf + '_bpf_env_out_' + str(i) + '.png'
            if save_png:
        	    plt.savefig(pathf)
        	    print ( pathf)
            else:
                plt.show()
                
                if 1: # set 1, save data as npz int the directory to use Noise-Spectrum-Estimation
                    pathf_npz= self.base_pathf + '_bpf_env_out_' + str(i) + '.npz'
                    np.savez(pathf_npz, bdata=d0, bfreq=f0, braw=self.fig[ :, p312[0] + i], \
                    bmax=np.amax( self.fig), bstep=(p312[0]+i))
                    print ('save a npz ', pathf_npz)
        
        
        # convert single Gray scale to RGB gray
        rgb_fig= np.zeros( (in_fig.shape[0],in_fig.shape[1], 3) )
        
        for i in range ( in_fig.shape[0] ):
            for j in range ( in_fig.shape[1] ):
                if in_fig[i,j] >= red_mark_threshold:  # RED mark if over red_mark_threshold
                    rgb_fig[i,j,0]= 1
                else:
                    rgb_fig[i,j,:]=  255 - in_fig[i,j]
        
        return rgb_fig


    def make_figure(self, sp=-1, ep=-1, pathf=None, save_png=False):
        # draw waveform
        if sp < 0 or ep < 0:
            sp=0
            ep=len(self.yg)
        
        delta_x= ep - sp
        delta_x_t= (delta_x / self.sr) * 1000.0 * 10
        self.fig= np.zeros( (self.num_bpf_bank, delta_x))
        self.wav_raw= np.zeros( (self.num_bpf_bank, delta_x))
        for i in range( self.num_bpf_bank):
            self.fig[i]= self.BPF_bank_wav_envlop[i][sp:ep].copy()
            self.wav_raw[i]= self.BPF_bank_wav[i][sp:ep].copy()
        
        # Normalize to [0, 255]
        if 1:  # linear
            f=  self.fig / np.amax( self.fig)
            self.fig_unit = np.uint8(np.around(f*255))
            
        else:  # log10
            f = np.log10( self.fig )
            #f=  f / np.amax( f )
            self.fig_unit = f #np.uint8(np.around(f*255))
        
        if pathf is not None:
            self.fig_rgb= self.gray2rgb(self.fig_unit, self.wav_raw, red_mark_threshold=999, show_number=0)
        else:  
            # if specify frequency range to detect max, only when frame number is specified
            if 0:
                print('warning: use specified frequency range to detect max.')
                fr0=[2000, 2600]
            else:
                fr0=None
            self.fig_rgb= self.gray2rgb(self.fig_unit, self.wav_raw, freq_rang=fr0, save_png=save_png)  
        
        
        spt= self.NSHIFT * self.frame_num / self.sr * 1000
        ept= (self.NSHIFT * self.frame_num  + delta_x) / self.sr * 1000
        
        fig, ax = plt.subplots()
        plt.title( 'frame no. ' + str(self.frame_num) +' : ' + str(spt) + ' to ' + str(ept) + ' msec'  )
        plt.imshow( self.fig_rgb, origin='lower', extent=[0, delta_x, self.Start_Freq, self.End_Freq ],\
        aspect=0.5) #1.0)
        
        #
        if pathf is not None:
        	plt.savefig(pathf)
        elif save_png:
            pathf2 = self.base_pathf + '_.png'
            plt.savefig(pathf2)
            print (pathf2)
        else:
            plt.show()
        
        
def get_path_name( dir0, path0, number0):
	# return file path name
    # make dir if the directory dir0 is not exist
    if not os.path.isdir( dir0 ):
        os.mkdir( dir0 )
    # get path0 basename without ext
    f, _= os.path.splitext( os.path.basename(path0))
    
    return dir0 + '/' + f + '_' + str(number0) + '.png'

def get_path_name2( dir0, path0, frame_num):
	# return file path name
    # make dir if the directory dir0 is not exist
    if not os.path.isdir( dir0 ):
        os.mkdir( dir0 )
    # get path0 basename without ext
    f, _= os.path.splitext( os.path.basename(path0))
    
    return dir0 + '/' + f + '_'  + str(frame_num)

if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser(description='BPF bank analysis')
    parser.add_argument('--wav_file', '-w', default='wav/i_1-16k.wav', help='wav-file-name(mono,16bit)')
    parser.add_argument('--frame', '-f', type=int, default=-1, help='specify the frame number, set negative value if ignore')
    parser.add_argument('--result_dir', '-r', default='result', help='specify directory to save result image')
    parser.add_argument('--en', action='store_true', help='save result image')
    parser.add_argument('--leak', action='store_true', help='see spectrum leak caused by the envelop') # not available with --en
    
    args = parser.parse_args()
    
    # instance
    sub1= Class_sub2( args.wav_file, frame_num=args.frame, result_dir=args.result_dir , save_png=args.en,\
    spectrum_leak= args.leak)
    

