# -*- coding: utf-8 -*-
"""
Created on Thu Jul 07 15:29:44 2016

@author: OPdaTe
"""
import numpy as np
from matplotlib.cm import jet
import matplotlib.pyplot as plt
import sys


class SCA_max:
    
    def __init__(self,surfaceinstance):
        self.array = surfaceinstance.array
        self.h_components = None
        self.h_components_sep = None
        self.h_peaks_num = None
        self.h_spacing_arr = None
        self.h_spacing_arr_sep = None
        self.h_results = None
        self.v_components = None      
        self.v_components_sep = None
        self.v_peaks_num = None
        self.v_spacing_arr = None
        self.v_spacing_arr_sep = None
        self.v_results = None
        self.full_peaks_component = None 
        self.max_comp = 7
        self.stage_step = surfaceinstance.parameters.stage_step

    def __progress(self,count, total, suffix=''):
        bar_len = 60
        filled_len = int(round(bar_len * count / float(total)))
    
        percents = round(100.0 * count / float(total), 1)
        bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
        sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
        sys.stdout.flush()  # As suggested by Rom Ruben   

    def find_max_peaks_horb(self,arr,transform=False):
        final_array=[]
        arr = np.ma.masked_invalid(arr)
        if transform:
            arr = arr.T
        for row in arr:
            row = np.insert(row,0,np.min(arr))
            row = np.append(row,np.min(arr))
            crow = [v for v in row if v==v]
            newrow=[]
            for index,hight in enumerate(crow[1:-1]):
                if hight>crow[index] and hight>crow[index+2]:
                    newrow.append(hight)
                else:
                    newrow.append(np.nan)
    
            final_array.append(newrow)
        #in case we have some nan in the original surface we wanto to get an array
        #with the original dimension
        if np.isnan(arr.data).any():
            adj_array=[]
        #What we have to do is to put all the put the nan of the original arry in 
        #in their original place and check if we have to mantain the peak or put 
        #another nan instead
            for row,newrow in zip(arr.data,final_array):
                adj_row=[]
                counter=0
                for value in row:
                    if np.isnan(value):
                        adj_row.append(np.nan)
    
                    else:
                        if np.isnan(newrow[counter]):
                            adj_row.append(np.nan)
                            counter+=1
                        else:
    #                        print 'Added %s' %(newrow[counter])
    #                        print counter
                            adj_row.append(newrow[counter])
                            counter+=1
                            
                #print adj_row   
                adj_array.append(adj_row)
            if transform: #Transform Back to the original shape
                adj_array = np.array(adj_array).T
            else:
                adj_array = np.array(adj_array)
                
            return np.array(adj_array)
        if transform:#Transform Back to the original shape
                final_array = np.array(final_array).T
        else:
                final_array = np.array(final_array)
                
        return final_array
    
    
    def Iter_peak(self,transform=False):
        surf = self.array
        peaks_number=[]
        surfaces=[]
        components=0
        kind = 'H'
        if transform:
            kind = 'V'
    
        while np.count_nonzero(~np.isnan(surf))>0 and components <self.max_comp:
            self.__progress(components, self.max_comp,
                            suffix='C: %s %s' %(components,kind))
            surf = self.find_max_peaks_horb(surf,transform)       
            peaks_number.append(np.count_nonzero(~np.isnan(surf)))
            surf= np.ma.masked_invalid(surf)
            surfaces.append(surf)
            components+=1
        masked_surfaces=[]       
        for index, surf  in enumerate(surfaces[:-1]):
            jsurf= surf.copy()
            xsurf = np.ma.mask_or(jsurf.mask,~surfaces[index+1].mask)
            jsurf.mask=xsurf
            masked_surfaces.append(jsurf)
            
        #retrieve results
        hist=[]
        percent=[]
        comulative=[]
        for index,value in enumerate(peaks_number[:-1]):
            res=value-peaks_number[(index+1)]
            hist.append(res)
            percent.append(float(res)/peaks_number[0]*100)
            comulative.append(sum(hist))
            
        if transform:
    
            self.v_components = surfaces      
            self.v_components_sep = masked_surfaces
            self.v_peaks_num = peaks_number
            self.v_results = [hist,percent,comulative]
 
        else:
            self.h_components = surfaces
            self.h_components_sep = masked_surfaces
            self.h_peaks_num = peaks_number
            self.h_results = [hist,percent,comulative]
             
    
            
        
        
    def plot_results(self):
        if self.v_results == None and self.h_results == None:
            print('Analysis not performed!')
            return
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        h_reslen0=0
        v_reslen0=0
        #it plots the results 
        #[0] histogram
        #[1] percent
        #[2] comulative 
        w = 0.3
        if self.v_results != None:
            if self.v_results[0] != []:
                ax3.bar(np.arange(1,len(self.v_results[0])+1),
                        self.v_results[0],width = w,
                         align='center',color='red',label='Vertical component')
                v_reslen0 = len (self.v_results[0])
        if self.h_results != None:
            if self.h_results[0] != []:
                ax3.bar(np.arange(1+w,len(self.h_results[0])+1+w),
                        self.h_results[0],width = w,
                         align='center',color='blue',label='Horizontal component')
                h_reslen0 = len (self.h_results[0])
                print('nofndofi')
        maxlen = max([h_reslen0,v_reslen0])  
        ax3.set_xticks(list(range(1,maxlen+1)))
        ax3.set_title('Number of peaks')
        ax3.set_xlabel(' Component ')
        ax3.legend()
        
        
        h_reslen1=0
        v_reslen1=0
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        if self.v_results != None:
            if self.v_results[1] != []:
                v_reslen1 = len (self.v_results[0])
                ax4.plot(np.arange(1,v_reslen1+1),self.v_results[1],'-x',label='Vertical component')
        if self.h_results != None:
            if self.h_results[1] != []:
                h_reslen1 = len (self.h_results[0])
                ax4.plot(np.arange(1,h_reslen1+1),self.h_results[1],'-x',label='Horizontal component')
            
        
        maxlen = max([h_reslen1,v_reslen1])
        ax4.set_xticks(list(range(1,maxlen+1)))
        ax4.set_title('Pertange of peaks for every component')
        ax4.set_ylabel(' Percentage of peaks (%)')
        ax4.set_xlabel(' Component ')
        ax4.legend()        
        ax4.set_xlim(0,maxlen+1)
        ax4.set_ylim(0,100)          
        
        h_reslen2=0
        v_reslen2=0
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111)
        if self.v_results != None:
            if self.v_results[2] != []:
                v_reslen2 = len (self.v_results[0])
                ax5.plot(np.arange(1,v_reslen2+1),self.v_results[2],label='Vertical component')
            
        if self.h_results != None:
            if self.h_results[2] != []:
                h_reslen2 = len (self.h_results[0])
                ax5.plot(np.arange(1,h_reslen2+1),self.h_results[2],label='Horizontal component')
                
        
        maxlen = max([h_reslen2,v_reslen2])
        ax5.set_xticks(list(range(1,maxlen+1)))
        ax5.set_xlim(0,maxlen+1)
        ax5.set_title('Comulative curve')
        ax5.set_ylabel(' Number of peaks')
        ax5.set_xlabel(' Component ')
        ax5.legend(loc=4)   
        
    def plot_surfaces(self,surfaces):
        figs={}
        axs={}
        cmap = jet
        cmap.set_bad('w',1.) #set bad values
        for idx,plot in enumerate(surfaces):
            figs[idx]=plt.figure()
            axs[idx]=figs[idx].add_subplot(111)
            axs[idx].pcolormesh(plot,cmap=cmap)
            plt.title('%s$^{th}$ Component' %(idx+1))
            
    def get_spacing_arrays(self,surfaces,spatial=False,transpose=False):
        #spacing array is useful to see if the space between the peaks is 
    #the same across the sample. 
        distance_arrays=[]
        for surf in surfaces:
            if transpose:
                surf = surf.T
            distance_array=[]
            for row in surf:
                indexes = np.ma.nonzero(row) #we get the indexes
                distance_array.append(np.diff(indexes)) #we caluclate the dist.
            distance_arrays.append(distance_array)
        
        adj_arrays=[]
        if spatial:
            for distarr,surf in zip(distance_arrays,surfaces):
                adj_array=[]
            #What we have to do is to put all the nan of the original arry in 
            #in their original place and check if we have to mantain the peak or put 
            #another nan instead
                for row,newrow in zip(surf,distarr):
                    adj_row=[]
                    counter=0
                    for value in row:
                        if np.isnan(value):
                            adj_row.append(np.nan)
                            print('app nan')
        
                        else:
                            if counter<newrow[0].size:
                                
                                if np.isnan(newrow[0][counter]):
                                    adj_row.append(np.nan)
                                    counter+=1
                                else:
        #                           print 'Added %s' %(newrow[counter])
        #                           print counter
                                    adj_row.append(newrow[0][counter])
                                    print(newrow[0][counter])
                                    counter+=1
                            else:
                                adj_row.append(np.nan)
                                
                    #print adj_row   
                        adj_array.append(adj_row)
                    adj_arrays.append(np.array(adj_array))
            if transpose:
                #Arrays with irregular shape are not easy to traspose
                Trasp = list(map(list, list(*distance_array)))
                distance_array = [x for x in Trasp if x is not None]
                return np.array(adj_arrays).T,distance_arrays
            else:
                return np.array(adj_arrays),distance_arrays
        
        else:
            
            if transpose:
                Trasp = list(map(list, list(*distance_array)))
                distance_array = [x for x in Trasp if x is not None]
                return distance_arrays
            else:
                return distance_arrays
     
    
    def get_spacing_results(self,spacing_arrays,kind=''):
        figs={}
        axs={}
        means=[]
        stdevs=[]
        #from matplotlib import rc
        #rc('text', usetex=True)    
        for idx, spc_arr in enumerate(spacing_arrays):
            idxt = str(idx)+'hist'
            res = np.hstack(spc_arr)[0]*self.stage_step
            if res.size != 0:
                figs[idxt]=plt.figure()
                bins = np.linspace(min(res),max(res),75)
                axs[idxt]=figs[idxt].add_subplot(111)
                axs[idxt].hist(res,bins=bins)
                axs[idxt].set_title(' %s $^{th}$ %s Component' %(idx+1,kind))
                axs[idxt].set_xlabel('Distance between peaks (mm)')
                mean = np.mean(res)
                stdev = np.std(res)
                means.append(mean)
                stdevs.append(stdev)
    #            x = axs[idx].get_xlim()[1]*0.8
    #            y = axs[idx].get_ylim()[1]*0.8
                table = 'Mean: %0.2f \nStd.Dev.: %0.2f' %(mean,stdev)
                axs[idxt].text(0.8,0.8,table,size=12,
                            transform = axs[idxt].transAxes,
                            bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

    def locate_full_peaks(self):
        self.full_peaks_component =(self.v_components[0] + self.h_components[0])/2
        
   
    def perform_component_analysis(self):
        self.Iter_peak(transform=False)
        print('Trasposing the array...')
        self.Iter_peak(transform=True)
        print('Locating full peaks..')
        self.locate_full_peaks()
        #self.get_spacing_arrays(spatial=False)
    
    def compute_spacing_array(self,sep = False):
        if sep:
            self.h_spacing_arr_sep = self.get_spacing_arrays(self.h_components_sep)
            self.v_spacing_arr_sep = self.get_spacing_arrays(self.v_components_sep,transpose=True)
        else:
            self.h_spacing_arr = self.get_spacing_arrays(self.h_components)
            self.v_spacing_arr = self.get_spacing_arrays(self.v_components,transpose=True)
     
    def plot_spacing_results(self, sep =False):
        if sep:            
            self.get_spacing_results(self.v_spacing_arr_sep,kind='vertical')       
            self.get_spacing_results(self.h_spacing_arr_sep,kind='horizontal')
            
        else:
             self.get_spacing_results(self.v_spacing_arr,kind='vertical')       
             self.get_spacing_results(self.h_spacing_arr,kind='horizontal')  
            
   
   
class SCA_min(SCA_max):
    '''
    SCA_min is the same of SCA_max we overwrite find_max_peaks_horb so that it 
    takes the min values instead of the max
    '''
        
    
    def find_max_peaks_horb(self,arr,transform=False):
        final_array=[]
        arr = np.ma.masked_invalid(arr)
        if transform:
            arr = arr.T
        for row in arr:
            row = np.insert(row,0,np.min(arr))
            row = np.append(row,np.min(arr))
            crow = [v for v in row if v==v]
            newrow=[]
            for index,hight in enumerate(crow[1:-1]):
                #HERE IS THE MAIN CONDITION
                if hight<crow[index] and hight<crow[index+2]:
                    newrow.append(hight)
                else:
                    newrow.append(np.nan)
    
            final_array.append(newrow)
        #in case we have some nan in the original surface we wanto to get an array
        #with the original dimension
        if np.isnan(arr.data).any():
            adj_array=[]
        #What we have to do is to put all the put the nan of the original arry in 
        #in their original place and check if we have to mantain the peak or put 
        #another nan instead
            for row,newrow in zip(arr.data,final_array):
                adj_row=[]
                counter=0
                for value in row:
                    if np.isnan(value):
                        adj_row.append(np.nan)
    
                    else:
                        if np.isnan(newrow[counter]):
                            adj_row.append(np.nan)
                            counter+=1
                        else:
    #                        print 'Added %s' %(newrow[counter])
    #                        print counter
                            adj_row.append(newrow[counter])
                            counter+=1
                            
                #print adj_row   
                adj_array.append(adj_row)
            if transform: #Transform Back to the original shape
                adj_array = np.array(adj_array).T
            else:
                adj_array = np.array(adj_array)
                
            return np.array(adj_array)
        if transform:#Transform Back to the original shape
                final_array = np.array(final_array).T
        else:
                final_array = np.array(final_array)
                
        return final_array
    

class SCA:
    def __init__(self, surfobj):
        self.smin = SCA_min(surfobj)
        self.smax = SCA_max(surfobj)
        self.h_components = []
        self.v_components = []
        self.noise_treshold = None
    
    def run(self):
        print("##PERFORMING ANALYSIS ON MAXIMUMS... ")
        self.smax.perform_component_analysis()
        print("##PERFORMING ANALYSIS ON MINIMUMS... ")
        self.smin.perform_component_analysis()
        hmin = np.min((len(self.smin.h_components),len(self.smax.h_components)))
        vmin = np.min((len(self.smin.v_components),len(self.smax.v_components)))
        print("##COMBINIG THE RESULTS...")
        for i in range(hmin):
            maxmask = self.smax.h_components[i].mask
            newmask = self.smin.h_components[i].mask & maxmask
            newarr = self.smin.h_components[i].copy()
            newarr[~maxmask] = self.smax.h_components[i][~maxmask]
            newarr.mask = newmask
            self.h_components.append(newarr)
        for i in range(vmin):
            maxmask = self.smax.v_components[i].mask
            newmask = self.smin.v_components[i].mask & maxmask
            newarr = self.smin.v_components[i].copy()
            newarr[~maxmask] = self.smax.v_components[i][~maxmask]
            newarr.mask = newmask
            self.v_components.append(newarr)
    
    def plot_surfaces(self,surfaces):
        figs={}
        axs={}
        cmap = jet
        cmap.set_bad('w',1.) #set bad values
        for idx,plot in enumerate(surfaces):
            figs[idx]=plt.figure()
            axs[idx]=figs[idx].add_subplot(111)
            axs[idx].pcolormesh(plot,cmap=cmap)
            plt.title('%s$^{th}$ Component' %(idx+1))

    def viewer(self):
            # Based on PythonQwt https://pypi.python.org/pypi/guiqwt
            print('NOT IMPLEMENTED YET!')
            sx = 0
            sy = 0
    
            from guiqwt.plot import ImageDialog
            from guiqwt.builder import make
    
            from guiqwt.tools import (
                RectangleTool,
                EllipseTool,
                HRangeTool,
                PlaceAxesTool,
                MultiLineTool,
                FreeFormTool,
                SegmentTool,
                CircleTool,
                AnnotatedRectangleTool,
                AnnotatedEllipseTool,
                AnnotatedSegmentTool,
                AnnotatedCircleTool,
                LabelTool,
                AnnotatedPointTool,
                VCursorTool,
                HCursorTool,
                XCursorTool,
                ObliqueRectangleTool,
                AnnotatedObliqueRectangleTool)
    
            from guiqwt.styles import style_generator, update_style_attr
    
            STYLE = style_generator()
    
            def customize_shape(shape):
                param = shape.shapeparam
                style = next(STYLE)
                update_style_attr(style, param)
                param.update_shape(shape)
                shape.plot().replot()
    
            def create_window():
                win = ImageDialog(
                    edit=True,
                    toolbar=True,
                    wintitle=self.name,
                    options=dict(
                        show_contrast=True,
                        show_xsection=True,
                        show_ysection=True,
                        xlabel='stage res.= %s mm x' %
                        (self.parameters.stage_step),
                        xunit='mm',
                        ylabel="y",
                        yunit='mm',
                        zlabel='Distance',
                        zunit='mm',
                        show_itemlist=True))
                win.resize(800, 600)
                for toolklass in (
                        LabelTool,
                        HRangeTool,
                        VCursorTool,
                        HCursorTool,
                        XCursorTool,
                        SegmentTool,
                        RectangleTool,
                        ObliqueRectangleTool,
                        CircleTool,
                        EllipseTool,
                        MultiLineTool,
                        FreeFormTool,
                        PlaceAxesTool,
                        AnnotatedRectangleTool,
                        AnnotatedObliqueRectangleTool,
                        AnnotatedCircleTool,
                        AnnotatedEllipseTool,
                        AnnotatedSegmentTool,
                        AnnotatedPointTool):
                    win.add_tool(toolklass)
                return win
    
            import guidata
            _app = guidata.qapplication()
            # --
            win = create_window()
            
            #plotting max h
            for i in self.smax.h_components:
                    
                array = make.image(
                    self.array,
                    title="Data",
                    interpolation="nearest")
        
                plot.add_item(array, z=3)
                win.exec_()
                plt.close('all')
                return win





if __name__ == '__main__':
    
    def sinusoidal_surf(finalx=10,res=0.1):
        Y=[np.sin(x**2)*5 for x in np.arange(0,finalx,res)]
        Z=[]
        for i in range(int(finalx/res)/2):
            Z.append(Y)
        return np.array(Z,dtype=np.float32)
    
    class param:
        def __init__(self):
            self.stage_step = 1
    
    class Test:
        #for compatibility with micro PRO
        def __init__(self,array):
            self.array = array
            self.parameters = param()
            
    Testsurf = Test(sinusoidal_surf())
    test = SCA_max(Testsurf)    