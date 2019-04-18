#coding:utf-8

# A trial estimation of noise source spectrum condition by fitting Gaussian distribution from BPF bank output
#
# BPF出力をガウス分布に当てはめることによりノイズのスペクトルの状態を予想する

import argparse
import glob
import numpy as np
from scipy import signal
from scipy import stats
from scipy.stats import norm
from scipy import optimize
from matplotlib import pyplot as plt


# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  matplotlib  2.1.1
#  scipy 1.0.0

class Class_estimate_nss1(object):
    def __init__(self, path0, top_level=255, show0=False ):
        # load result data (BPF output)
        self.path0= path0
        xdata= np.load(self.path0)
        self.d0= xdata['bdata']
        self.f0= xdata['bfreq']
        self.fig_raw= xdata['braw']
        self.fig_max= xdata['bmax']
        self.d1= (self.fig_raw / self.fig_max) * top_level
        self.p312i= xdata['bstep']
        print ('load result data from ', self.path0)
        
        # get  peak candidate
        self.freq_min=None
        self.f_result, self.i_result= self.pk_detect( self.d1, self.f0)
        # get  drop-peak candidate
        self.f_result2, self.i_result2= self.pk_detect( top_level -self.d1 , self.f0)
        
        # get fitting points candidate
        self.v_result3, self.i_result3, index_maxpk, f_peak= self.get_maxpk_around(self.d1, self.f0, \
        self.i_result, self.i_result2, db_min=6)
        
        # get suitable parameter by downhill simplex algorithm.
        X_init=[f_peak, index_maxpk, 30.0]
        X_value=self.i_result3
        Y_target=self.v_result3
        
        res0= self.minimize_cost(X_init, X_value, Y_target)  # return estimated parameter
        
        self.A=res0[0]
        self.loc=res0[1]
        self.scale=res0[2]
        self.delta_freq= self.f0[1] - self.f0[0]
        self.freq_center= self.f0[0] + self.delta_freq * self.loc
        self.freq_stddev= self.delta_freq * self.scale
        print ('')
        print ('* estimation result *')
        print ('freq center [Hz] ', self.freq_center)
        print ('standard deviation ', self.freq_stddev)
        
        # get strength
        self.i_result_c_maxsort, self.strengths= self.get_strength()
        
        # show estimation result
        x0= np.arange( len(self.d1))
        y=self.gauss( x0, self.A, self.loc, self.scale)
        
        if show0:
            self.plot_figure0( self.d1, self.f0, self.i_result, self.i_result2, self.i_result3, y, \
             self.i_result_c_maxsort, self.strengths)
        
        
    def gauss(self, x, A, loc, scale):
        # The location (loc) keyword specifies the mean. 
        # The scale (scale) keyword specifies the standard deviation.
        # norm.pdf(x) = exp(-x**2/2)/sqrt(2*pi)
        # norm.pdf(x, loc, scale)= norm.pdf(y) / scale with y = (x - loc) / scale.
        return  (A * scale * np.sqrt(2 * np.pi))  * stats.norm.pdf(x, loc, scale)

    def minimize_cost(self, X_init, X_value, Y_target, disp=False):
        # try to minimize the function
        #   by "fmin" that is minimize the function using the downhill simplex algorithm.
        args1=(X_value, Y_target)
        res_brute = optimize.fmin( self.calc_cost, X_init, args=args1, full_output=True, disp=disp)
        print ( 'minimum ', res_brute[0] )  # minimum
        print ( 'function value ', res_brute[1] )  # function value at minimum
        if res_brute[4] != 0:  # warnflag
            print ('warnflag is not 0')
        return res_brute[0]

    def cost_0(self, y_target, y):
        # lower cost function
        return abs(y_target - y).mean()

    def calc_cost(self, X , x_value, y_target):
        # get mean of difference between target and new computed ones
        y= self.gauss( x_value, X[0], X[1], X[2])
        cost0= self.cost_0(y_target, y)
        return cost0


    def get_strength(self,):
        # ピーク周波数の候補ポイントにおいて
        # 推定されたガウス分布の値より、何ｄＢ大きいか（=strength, 強さ）を求める
        #
        #   出力：ピーク周波数の候補ポイントで大きいもの順のインデックス
        #         推定されたガウス分布から差　ｄＢ大きいか（=strength, 強さ）
        
        i_result=np.array( self.i_result)
        i_result_c= i_result[ (i_result > self.i_result3[0]) & (i_result < self.i_result3[-1]) ]
        index_max_sort= np.argsort( self.d1[i_result_c])[-1::-1]
        i_result_c_maxsort= i_result_c[ index_max_sort]  # ピーク周波数の候補ポイントで大きいもの順index
        
        strengths= np.zeros( len(i_result_c_maxsort) )
        
        for i, x in enumerate(i_result_c_maxsort):
            y=self.gauss( x, self.A, self.loc, self.scale)
            strengths[i]= 20. * np.log10( self.d1[x] / y)  # strength unit dB
            if i==0:
                print ( 'frequency[Hz]    peak-value  strength[dB]' )
            #print ( self.f0[0] + self.delta_freq * x, self.d1[x], strengths[i] )
            print ( '{0}             {1:.1f}       {2:.1f} '.format(self.f0[0] + self.delta_freq * x,\
             self.d1[x], strengths[i]) )
        
        print ('')
        print ('')
        
        return i_result_c_maxsort, strengths

    def get_maxpk_around(self, d1, f0, i_result, i_result2, db_min=6):
        #　ピーク周波数の最大値を見つけ
        #  その周囲でdb_min dB以内のドロップピークを探す
        #　
        #  入力：スペクトル
        #       ピーク周波数リスト
        #       ドロップピーク周波数リスト
        #      （オプション）最大のピーク値の隣のドロップピークよりdb_min dB以内の箇所を高い候補とする
        #
        #   出力：候補のインデックス
        #         候補の値       
        #         最大のピークのインデックス
        #         最大のピーク値の隣のドロップピークの値
        
        v_result3_high=[]
        i_result3_high=[]
        v_result3_low=[]
        i_result3_low=[]
        
        # ピークを大きいもの順にならべ変える
        index_sort= np.argsort( d1[i_result ] )[-1::-1]
        index_max_sort= np.array(i_result)[index_sort]
        #print ( d1[ index_max_sort] )
        index_maxpk= index_max_sort[0]
        db_gain= np.power(10.0, -1.0 * db_min / 20.) # convert from dB to multiplier
        
        
        # 最大ピークより高い周波数の方を探す
        next_index_high=-1
        next_index_high_value=-1
        for i, v0 in enumerate( i_result2 ):
            if v0 > index_maxpk:
                next_index_high=i
                next_index_high_value= d1[v0]
                break
                
        # 最大ピークより低い周波数の方を探す
        next_index_low=-1
        next_index_low_value= -1
        for i, v0 in enumerate( i_result2 ):
            if v0 >= index_maxpk:
                break
            else:
                next_index_low=i
                next_index_low_value=d1[v0]
        
        # 最大のピーク値の隣のドロップピークよりdb_min dB
        f_peak= np.amax( (next_index_high_value, next_index_low_value) )
        f_peak_6dB= f_peak * db_gain
        
        if next_index_high > 0: # 高い周波数にドロップピークが存在するときは
            for i in range (next_index_high, len(i_result2)):
                if d1[ i_result2[i] ] >= f_peak_6dB:
                    v_result3_high.append( d1[ i_result2[i] ] )
                    i_result3_high.append( i_result2[i] )     
                else:  # 最大のピーク値からdb_min dB未満なので　より高い周波数の探索は終了する
                    break
        
        if next_index_low > 0: # 低い周波数にドロップピークが存在するときは
            for i in range (next_index_low, 0,-1):
                if d1[ i_result2[i] ] >= f_peak_6dB:
                    v_result3_low.append( d1[ i_result2[i] ] )
                    i_result3_low.append( i_result2[i] )     
                else:  # 最大のピーク値からdb_min dB未満なので　より低い周波数の探索は終了する
                    break
        
        # 周波数が低いもの順に並べ変える
        v_result3= v_result3_low[-1::-1] + v_result3_high
        i_result3= i_result3_low[-1::-1] + i_result3_high
        
        return list(v_result3), list(i_result3), index_maxpk, f_peak
        

    def pk_detect(self, d1, f0, db_min=3):
        #   山型（凸）のピークポイントを見つける
        #
        #   入力：スペクトル
        #         周波数リスト
        #         （オプション）周囲からdb_min dB以上高い場合をピークとみなす
        #
        #   出力：ピークのインデックス
        #         ピークの周波数
        is_find_first= False
        f_result=[]
        i_result=[]
        
        for i in range (1,len(d1)-1):
            if self.freq_min is not None and  f0[i] <= self.freq_min :
                continue
                
            if d1[i] > d1[i-1] and d1[i] > d1[i+1] :
                if not is_find_first :
                    f_result.append( f0[i])
                    i_result.append(i)
                    is_find_first =True
                else:
                    f_result.append( f0[i] )
                    i_result.append(i)
        
        
        # 周囲からdb_min dB以上高い場合をピークとみなす
        Q_list,_,_,_,_ =  self.bandwidth_detect( d1, f0, i_result, db_min=3)
        f_result2= np.array(f_result)[np.where( np.array(Q_list) > 0.0)[0]]
        i_result2= np.array(i_result)[np.where( np.array(Q_list) > 0.0)[0]]
        
        return list(f_result2), list(i_result2)  # 旧品との互換性を保つためリストに戻す


    def bandwidth_detect(self, d1, f0, peak_index_list, db_min=3):
        #   スペクトルのピークw0から
        #   3dB(db_minのデフォルト値)低下した周波数を求める
        #
        #   入力：対数スペクトル  * 高精度の求めるには　周波数分解能が高いことが求められる。
        #         周波数単位
        #         ピークのインデックス
        #
        #   出力：-3dB低下した周波数w1（ピークより低い周波数の方）
        #         -3dB低下した周波数w2（ピークより高い周波数の方）
        #          Q= w0/(w2-w1) 定義できないときは零が入る
        #          暫定的に　w1 又は w2のどちらかが分かっているときは、Q= w0/(2*(w0-wx))を入れておく
        Q_list=[]
        Low_freq_list=[]
        High_freq_list=[]
        Low_freq_index=[]
        High_freq_index=[]
        db_gain= np.power(10.0, -1.0 * db_min / 20.) # convert from dB to multiplier
        
        for ipk in (peak_index_list):
            low_index= 0
            high_index= 0
            f_peak= d1[int(ipk)]
            f_peak_3dB= d1[int(ipk)] * db_gain
            
            # ピークより低い周波数の方を探す
            for i in range (int(ipk),0,-1):
                if d1[i] <= f_peak_3dB:
                    low_index= i
                    break
                elif d1[i] > f_peak:
                    break
            
            # ピークより高い周波数の方を探す
            for i in range (int(ipk),len(d1)):
                if d1[i] <= f_peak_3dB:
                    high_index= i
                    break
                elif d1[i] > f_peak:
                    break
            
            #Q= w0/(w2-w1)の計算　 定義できないときは零が入る
            if low_index > 0 and high_index > 0:
                Q= 1.0 * int(ipk) / (high_index - low_index)
            elif low_index > 0:
                Q= 1.0 * int(ipk) / (2.0 * (int(ipk) - low_index))  # 暫定値を入れておく
            elif high_index > 0:
                Q= 1.0 * int(ipk) / (2.0 * ( high_index - int(ipk)))  # 暫定値を入れておく
            else:
                Q=0.0
            
            Q_list.append( Q )
            Low_freq_list.append( f0[low_index])
            High_freq_list.append( f0[high_index])
            Low_freq_index.append( low_index)
            High_freq_index.append( high_index)
            
            #print ( Q, f0[low_index], f0[high_index], f0[int(ipk)])
        
        return Q_list, Low_freq_list, High_freq_list, Low_freq_index, High_freq_index


    def plot_figure0(self, d0, f0, index_peak=None, index_drop_peak=None, index_fitting=None, y_gauss=None,\
     i_result_c_maxsort=None, strengths=None):
        # show spectrum (= BPF bank output )
        plt.figure()
        plt.title('BPF bank output: no.' + str( int(self.p312i) ))
        plt.xlabel('Hz')
        plt.plot( f0,  d0 )
        if index_peak is not None: 
            plt.plot( f0[index_peak], d0[ index_peak ] , "ro")
        else:
            #index_peak= signal.find_peaks_cwt( d0 , [ 2 ,3 ])
            pass
        if index_drop_peak is not None:
            plt.plot( f0[index_drop_peak], d0[ index_drop_peak ] , "bo")
            
        if index_fitting is not None:
            plt.plot( f0[index_fitting], d0[ index_fitting ] , "go" ) #color="c", marker="x")
        
        if y_gauss is not None:
            plt.plot( f0, y_gauss, color="r", ls="--")
        
        if i_result_c_maxsort is not None:
            for i,x in enumerate(i_result_c_maxsort):
                plt.text( f0[x],d0[x], format(strengths[i], '.1f') )
        
        plt.grid()
        plt.tight_layout()
        plt.show()
        
        
        
# misc helper function
def get_list(dir_in):
    # ファイル拡張子
    ext0='.npz'
    List0=glob.glob(dir_in + '/*' ,recursive=True)
    List1=[s for s in List0 if s.endswith( ext0 ) ]
    print ('number of files ', len(List1))
    return List1

if __name__ == '__main__':
    #
    parser = argparse.ArgumentParser(description='estimation noise spectrum condition')
    parser.add_argument('--npz_dir', '-r', default='../result_i', help='specify the directory of BPF bank data, *bpf_env_out*.npz ')
    args = parser.parse_args()
    
    # get result file list
    Flist1= get_list( args.npz_dir )
    
    # instance
    for path0 in Flist1:
        ana1=Class_estimate_nss1( path0 ,show0=True )
        