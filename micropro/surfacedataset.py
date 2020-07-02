from __future__ import print_function
import sys
if hasattr(__builtins__, 'raw_input'):
    input = raw_input
w = 'wb' # for csv module
a = 'ab'
if sys.version_info[0] >= 3:
    w = 'w'
    a = 'a'
# End python 2 backward compatibility
from .surface import ns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import jet, Greys
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from os import listdir, getcwd, path, makedirs,remove
from .config import myconfig
import time




__all__ = ['surfacedataset']

class surfacedataset:

    def __init__(self):
        self.surfaces = []
        self.name = getcwd().split('\\')[-1]
        self.log = []
        self.totaltime = None
        self.ID = None
        self.ID_b = None
        self.dataset_infos = None
        self.ROIs_dataset = None

    def save_log(self, name = None):
        if name == None:
            name = r'dataset_processing_log.txt'

        with open(name, 'w') as f:
            f.write('************************ \n')
            f.write('* Dataset: %s           *\n' % (self.name))
            f.write('************************ \n')
            for i in self.surfaces:
                f.write('> %s \n' %(i.name))
            f.write('************************ \n')
            f.write('* PROCESSING PERFORMED  *\n')
            f.write('************************ \n')
            for i in self.log:
                f.write('> %s \n' %(i))

    def load_data(self):
        pass

    def load_data_fromfolder(self,path = None):
        from os.path import isdir, join
        if path == None:
            path = getcwd()
            folders = listdir(path)
            for i in folders:
                if isdir(i) and i != 'Results' and i[0] != '.':
                    self.surfaces.append(ns(join(i,i)))
        else:
             folders = listdir(path)
             for i in folders:
                 j = join(path,i)
                 if isdir(j) and i != 'Results' and i[0] != '.':
                    finalpath = join(j,i)
                    print(finalpath)
                    self.surfaces.append(ns(finalpath))



    def load_mat_fromfolder(self,lateral_resolution,conversion_factor=1):
        import scipy.io
        matfiles = [i for i in listdir(getcwd()) if i[-4:] == '.mat']
        for i in matfiles:
            surf = ns(None,init=False)
            mat = scipy.io.loadmat(i)
            name = 'm'+i.split('.')[0]
            surf.array = np.array(mat[name])*conversion_factor
            surf.name = name
            surf.parameters.stage_step = lateral_resolution
            self.surfaces.append(surf)

    def load_path(self,path):
        self.surfaces.append(ns(path,init=False))

    def check_label_name(self):
        LabelsNotInName=0
        noIDs=0
        status=True
        datasetID = []
        for i in self.surfaces:
            if i.sample_infos != None:
                label = i.sample_infos.label
                name = i.sample_infos.name
                idx = i.header['JOB_ID']
                suff = ''
                datasetID.append(i.sample_infos.dataset)
                if  label not in i.name:
                    LabelsNotInName+=1
                    suff = "  WRONG!->"
                print("%s%s: %s  (%s ID:%s) " %(suff,i.name,label,name,idx))

            else:
                print(" %s  has no ID!" %(i.name))
                status = False
                noIDs+=1
        if datasetID != []:

            self.ID = list(set(datasetID))
            print("Found %s dataset:" %(len(self.ID)))
            try:
                import LabDBwrapper as Lab
                p = myconfig['databasepath']
                if p is not None and p != '':
                    self.dataset_infos = Lab.dataset(self.ID[0],p)
            except ImportError:
                "No LabWrapper found"

            print(self.ID)
        if not status or LabelsNotInName:
            print("******WARNING!!!!********")
            print("Samples without label: %s" %(noIDs))
            print("Samples with label not in name: %s" %(LabelsNotInName))
            if input("Continue? y/n").upper() == "Y":
                status = True
        return status

    def get_surface_by_name(self,name):
        surfaces_index = []
        for index,surf in enumerate(self.surfaces):
            if surf.name == name:
                surfaces_index.append(index)
        return surfaces_index

    def get_surface_by_ID(self,ID):
        surfaces_index = []
        for index,surf in enumerate(self.surfaces):
            if surf.sample_infos != None:
                if surf.sample_infos.id == ID:
                    surfaces_index.append(index)
            else:
                print("NO samples infos!")
        return surfaces_index

    def get_surface_by_label(self,label):
        surfaces_index = []
        for index,surf in enumerate(self.surfaces):
            if surf.sample_infos != None:
                if surf.sample_infos.label == label:
                    surfaces_index.append(index)
            else:
                print("NO samples infos!")
        return surfaces_index

    def dict_surfaces_by_name(self):
        """
        Return a dictionary with all the surfaces
        """
        dict_surfaces = {}
        for i in self.surfaces:
            dict_surfaces[i] = i.name
        return dict_surfaces

    def check_status(self,ignoreQC=False):
        self.check_label_name()
        self.get_QualityControl()
        done = 0
        if self.dataset_infos != None:
            for i in self.dataset_infos.samples:
                 print("Sample: %s :" %(i[1]))
                 indexs = self.get_surface_by_ID(i[0])
                 if indexs == []:
                     print(' *** MISSING ***')
                 else:
                     res = False
                     if ignoreQC:#Sample is counted as done even if QC fail
                         res=True
                     for j in indexs:
                         space = '          '+' '*len(i[1])
                         name = self.surfaces[j].name
                         QC = self.surfaces[j].logQC['result']
                         print("%s %s QC: %s" %(space,name,QC))
                         if QC == 'Passed':
                             res = True
                     if res:
                         done += 1

            self.time_evaluation()
            status = done /float(self.dataset_infos.totnumber)
            print("STATUS:{0:.0f}%".format(status*100))
            timeleft = (self.totaltime/(int(status*100)))*int((100 -status*100))
            print("TIME LEFT: %s" %(timeleft))

        else:
            print("No dataset infos!")

    def fliprows(self):
        for i in self.surfaces:
            print(i.name)
            i.fliprows()
            print("------------------")
        self.log.append("Flipped rows")

    def X_VelTest(self):
        pass
    def combine_usingMASK(self,indices = [],deletesurf= False,addsurf=True):
        surf2combine = []
        arr2combine = []
        laserpower = []
        for i in indices:
            surf2combine.append(self.surfaces[i])
            arr2combine.append(self.surfaces[i].array)
            laserpower.append( self.surfaces[0].parameters.laserpower)
        arr = np.ma.array(arr2combine)
        res = np.ma.mean(arr, axis = 0)
        import copy
        newsurf = copy.deepcopy(self.surfaces[0])
        newsurf.array = res
        newsurf.parameters.laserpower = set(laserpower)
        if deletesurf:
            for i in sorted(indices,reverse = True):
                del self.surfaces[i]
        if addsurf:
            self.surfaces.append(newsurf)
        repscore = np.ma.std(arr, axis = 0)
        return newsurf,repscore


    def combine_usingSNR(self,i=0,j=1,getSNR=True,threshold=500):
        surfa = self.surfaces[i]
        surfb = self.surfaces[j]
        laserpowerdiff = surfa.parameters.laserpower-surfb.parameters.laserpower
        print('Laser power difference: %s'%(laserpowerdiff))
        if surfa.parameters.laserpower < surfb.parameters.laserpower:
            bestsurf = surfa
            worsesurf = surfb

        else:
            bestsurf = surfb
            worsesurf = surfa

        if surfa.array.shape == surfb.array.shape:
            SNRbestsurf = bestsurf.SNR(matrix=True,plot=False)[1]
            if threshold != None:
                #if we set a threshold all the value grater than the threshold
            #will be bring to 1000 so that the bestmatrix will be kept
                whereisgreaterofthreshold = SNRbestsurf > threshold
                SNRbestsurf[whereisgreaterofthreshold]=1000

            SNRworsesurf = worsesurf.SNR(matrix=True,plot=False)[1]
            SNRaisbest = SNRbestsurf > SNRworsesurf#map where SNRa is best
            newarray = np.zeros(surfa.array.shape)
            newarray[SNRaisbest] = bestsurf.array[SNRaisbest]
            newarray[~SNRaisbest] = worsesurf.array[~SNRaisbest]
            snrarray=None
            if getSNR:
                snrarray = np.zeros(bestsurf.array.shape)
                snrarray[SNRaisbest] = SNRbestsurf[SNRaisbest]
                snrarray[~SNRaisbest] = SNRworsesurf[~SNRaisbest]
            import copy
            newsurf = copy.deepcopy(bestsurf)
            newsurf.array = newarray

            return newarray,snrarray,newsurf
        else:
            print('Shapes do not match!')

    def compute_Roughenss_ROI(self):
        from scipy.stats.mstats import skew, kurtosis
        import csv
        fields = [
            'Name',
            'Sq1',
            'Sq2',
            'Sq3',
            'Mean',
            'Skw1',
            'Skw2',
            'Skw3',
            'Mean',
            'Kurt1',
            'Kurt2',
            'Kurt3',
            'Mean']
        with open(r'Sq_comparison.csv', a) as f:
            writer = csv.writer(f)
            writer.writerow(fields)

            for surf in self.surfaces:
                surf.mfilter()
                Sq = []
                Sk = []
                Kur = []
                dROI1 = surf.array[400:550, 500:650]
                dROI1 = surf.subtractplane(array=dROI1)
                SqROI1 = round(np.std(dROI1) * 1000, 2)
                SkROI1 = round(skew(dROI1.flatten()), 2)
                KuROI1 = round(kurtosis(dROI1.flatten()), 2)
                Sq.append(SqROI1)
                Sk.append(SkROI1)
                Kur.append(KuROI1)

                dROI2 = surf.array[400:550, 750:900]
                dROI2 = surf.subtractplane(array=dROI2)
                SqROI2 = round(np.std(dROI2) * 1000, 2)
                SkROI2 = round(skew(dROI2.flatten()), 2)
                KuROI2 = round(kurtosis(dROI2.flatten()), 2)
                Sq.append(SqROI2)
                Sk.append(SkROI2)
                Kur.append(KuROI2)

                dROI3 = surf.array[400:550, 1250:1400]
                dROI3 = surf.subtractplane(array=dROI3)
                SqROI3 = round(np.std(dROI3) * 1000, 2)
                SkROI3 = round(skew(dROI3.flatten()), 2)
                KuROI3 = round(kurtosis(dROI3.flatten()), 2)
                Sq.append(SqROI3)
                Sk.append(SkROI3)
                Kur.append(KuROI3)

                average = np.mean(Sq)
                Sq.append(average)
                averages = np.mean(Sk)
                Sk.append(averages)
                averagek = np.mean(Kur)
                Kur.append(averagek)
                writer.writerow([surf.name] + Sq + Sk + Kur)

    def plot(self, vminx=None, vmaxx=None, unit='mm', data='', absolute=False,
             mode=False, plot=True,save=False,cartesian=True,
            dilution=False, title = True, cm = None,
            QC = False ,orientation = 'auto',show_ROIs=False):

        for i in self.surfaces:
            i.plot(vminx=vminx, vmaxx=vmaxx, unit=unit, data=data,
                   absolute=absolute,mode=mode, plot=plot,
                   save=save,cartesian=cartesian,dilution=dilution,rectangle=[],
                   title = title,
                   cm = cm,QC = QC ,orientation = orientation,show_ROIs=show_ROIs)

    def profile_plot(self, col=None,row=None, start =None, stop =None,vminx=None,
                     vmaxx=None, unit='mm', absolute=False,mode=False,
                     plot=True,save=False, title = True):
        '''
        plot profile
        '''
        datas =[]
        x_datas = []
        labels = []
        kind = ''
        if col != None and row == None:

            for i in reversed(self.surfaces):
                if col == 'mid':
                    col = i.array.shape[1]/2
                profile = i.array[:,col][start:stop]
                datas.append(profile)
                labels.append(i.name)
                x_datas.append(np.arange(len(profile))*i.parameters.stage_step)
                kind = 'vertical'

        elif row !=None and col == None:

             for i in reversed(self.surfaces):
                if row == 'mid':
                    row = i.array.shape[0]/2

                profile = i.array[row][start:stop]
                datas.append(profile)
                labels.append(i.name)
                x_datas.append(np.arange(len(profile))*i.parameters.stage_step)
                kind = 'horizontal'

        else:
            print("Select a row or a column")


        #This stretch colorbar values to min-max
        if mode == 'min-max':
            vminx = np.array(datas).min()
            vmaxx = np.array(datas).max()
        #You can also stre
        if mode == 1 or mode == 2 or mode == 3:
            stdar, meanar = np.std(datas.flatten()), np.mean(data.flatten())
            vminx = meanar - mode * stdar
            vmaxx = meanar + mode * stdar

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for profile,x_data,name in zip(datas,x_datas,labels):
            if unit == 'micron':
                profile = np.array(profile)*1000
            ax.plot(x_data,profile,label = name)

        ax.legend()
        ax.set_xlabel('(mm)')
        ax.set_ylabel('(%s)' %(unit))
        ax.legend()
        if title:
            ax.set_title('%s %s' % (kind,self.name))

        if save: self._save_plot(fig,'%s profile com %s.png' %(kind,self.name))
        if plot:
            plt.show()
        else:
            return plt



    def time_evaluation(self):
        import datetime
        t0 = datetime.timedelta()
        for i in  self.surfaces:
            duration = i.parameters.time['duration']
            h,m,s = list(map(int,duration.split(':')))
            t0 += datetime.timedelta(hours = h,minutes = m,seconds = s)
        print(t0)
        totsec = t0.total_seconds()
        h = totsec // 3600
        m = (totsec %3600) // 60
        sec = (totsec %3600) %60
        self.totaltime = t0
        return h,m,sec

    def save_metrology_complete(self,filename= None,ROI=None,plot=True,cutlabel=None,
                                groupbylabel=False,margin=None,save=False,
                                split=False,delimiter=None,uselabel=False,
                                centerROI = False, offset = False, dim = False,
                                subtractplane = False):
        import csv
        fields = list(self.surfaces[0].calculate_area_metrology(margin=margin).keys())
        labels = []
        values = []
        class_maker = {}
        i=1
        if filename == None:
            filename = r'Metrology_complete.csv'
        with open(filename, a) as f:
            writer = csv.writer(f)
            if split != False:
                name = ['Name','Type','Ref']
            else:
                if uselabel or cutlabel:
                    name = ['Name','Label']
                else:
                    name = ['Name']
            writer.writerow(name+fields)

            for surf in self.surfaces:
                #calculate the data using the methodo of ns class
                label = ' '
                if centerROI:
                    ycenter,xcenter = np.array(surf.array.shape)/2
                    y1 = offset + ycenter - dim
                    y2 = offset + ycenter + dim
                    x1 = offset + xcenter - dim
                    x2 = offset + xcenter + dim
                    ROI = [y1,y2,x1,x2]
                diz2 = surf.calculate_area_metrology(ROI = ROI,
                                                     margin = margin,
                                                     subtractplane = subtractplane)

                if (surf.sample_infos != None) and (cutlabel==None):
                    label = surf.sample_infos.label
                    labels.append(label)

                    if split == 'label':
                        writer.writerow([surf.name,label[0],label[1]]+list(diz2.values()))
                    elif split == 'filename':
                        if delimiter !=None:
                            dataclass = surf.name.split(delimiter)
                        else:
                            dataclass = surf.name

                        writer.writerow([surf.name,dataclass[-2],label[-1]]+list(diz2.values()))
                    else:
                        if uselabel:
                            writer.writerow([surf.name,label]+list(diz2.values()))
                        else:
                            writer.writerow([surf.name]+list(diz2.values()))

                    if groupbylabel:
                        if surf.sample_infos.label not in class_maker:
                            class_maker[surf.sample_infos.label] = i
                            i+=1
                else:
                    if cutlabel == None and split == False:
                        labels.append(surf.name)
                        writer.writerow([surf.name]+list(diz2.values()))

                    elif split == 'filename':
                        if delimiter !=None:
                            dataclass = surf.name.split(delimiter)
                        else:
                            dataclass = surf.name
                        writer.writerow([surf.name,dataclass[-2],label[-1]]+list(diz2.values()))

                    elif cutlabel != None:
                        label = surf.name[cutlabel[0]:cutlabel[1]]
                        labels.append(label)
                        writer.writerow([surf.name,label]+list(diz2.values()))

                values.append(list(diz2.values()))

        if plot:
            print("Plotting the data it will take some time plt.show() to show")
            figs={}
            axs={}
            values = np.array(values)
            for idx, field in enumerate(fields):
                figs[idx]=plt.figure()
                Y = values[:,idx]
                if groupbylabel:
                    X = [class_maker[j] for j in labels]
                else:
                    X = np.arange(1,len(labels)+1)
                xlim=max(X)+1
                axs[idx]=figs[idx].add_subplot(111)
                axs[idx].plot(X,Y,'x')

                axs[idx].hlines(np.mean(Y),0,xlim,label='Mean')
                axs[idx].set_xticks(X)
                axs[idx].set_title(field)
                axs[idx].set_ylabel('%s ' %(field))
                axs[idx].set_xlim(0,xlim)
                axs[idx].set_ylim(min(Y)*1.05,max(Y)*1.05)
                if groupbylabel:
                    axs[idx].set_xticklabels(np.unique(labels))

                else:
                    axs[idx].set_xticklabels(labels)
                axs[idx].grid(True)
                axs[idx].legend()
                if save:
                    self._save_plot(figs[idx],'%s Dataset ' %(field))

    def multiple_ROI_metrology(self,nROI = 3, offset = 100, dim = 200 ,
                               subtractplane = False):
        ty = (nROI -1)/2

        rang = list(range(-ty, ty+2))
        rang.remove(0)
        for i in rang:
            self.save_metrology_complete(centerROI = 1, offset = offset*i,
                                         dim = dim, plot=0,
                                         subtractplane = subtractplane )


    def workflow(self, flip_rows=True, checkmissing=True, mfilter=True,bilateral=False,
                 subtractplane = True, secmfilter = False, appyreport=False,
                 timingdiagram = False,saveres = False,
                 csvreport=True,cropconcentric=None,crop=None,
                 order=1,flipyax=None,metrology=3,checklabel=True,total='optimal'):
        import gc
        if checklabel:
            if self.check_label_name() == False:
                return False

        if flip_rows:
            print('FLIPPING ROWS OF THE SAMPLES:')
            self.fliprows()

        log = True
        for j, i in enumerate(self.surfaces):
            print('SURFACE %s of %s' % (j + 1, len(self.surfaces)))
            print('Name: %s' %(i.name))
            if checkmissing:
                print('Checking for missing values...')
                i.checkmissing(replacewith=0)
                if log: self.log.append('Checked missed value')

            if mfilter:
                print('Filtering bad values...')
                i.mfilter(total=total)
                if log: self.log.append('Bad value filter')

            if crop:
                print('Cropping concentric the matrix...')
                i.crop(crop[0],crop[1],crop[2],crop[3])
                if log: self.log.append('Cropped. Side: y1 %s y2 %s x1 %s x2 %s'
                %(crop[0],crop[1],crop[2],crop[3]))

            if cropconcentric !=None:
                print('Cropping concentric the matrix...')
                i.crop_concentric_ROI(cropconcentric)
                if log: self.log.append('Square Cropped. Side: %s mm' %(cropconcentric*2))

            if bilateral:
                print('Using bilateral filter ...')
                i.bilateralFilter()
                if log: self.log.append('Bilateral filter')

            if subtractplane:
                print('Subtracting plane...')
                i.subtractplane(order=order)
                if log: self.log.append('Plane subtraction')

            if secmfilter != False:
                print('Applying second filter...')
                i.mfilter('-%s' %(secmfilter),'+%s' %(secmfilter),snr=False)
                if log: self.log.append('Second bad value filter:-%s +%s' %(secmfilter))

            if timingdiagram:
                print('Timing diagram...')
                pltX = i.timing_diagram(plot=False)
                if saveres:
                    self._save_plot(pltX,'%s TD' %(i.name))
            if appyreport:
                print('Exporting report ...')
                i.exp_appyreport(flipyax=flipyax,metrology=metrology)
            gc.collect()
            gc.collect()
            log = False

            if saveres:
                pltFIG = i.plot(plot=False)
                pltSNR = i.SNR(plot=False)
                self._save_plot(pltFIG,'%s FIG' %(i.name))
                self._save_plot(pltSNR,'%s SNR' %(i.name))

        if csvreport:
            self.exp_csvReport()

    def viewer(self, use_local_ref=False):
        # Based on PythonQwt https://pypi.python.org/pypi/guiqwt
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
        plot = win.get_plot()

        if use_local_ref:
            for j, i in enumerate(self.surfaces):
                array = make.image(
                    i.array,
                    title=i.name,
                    xdata=[
                        0,
                        i.array.shape[1] *
                        i.parameters.stage_step],
                    ydata=[
                        0,
                        i.array.shape[0] *
                        i.parameters.stage_step],
                    interpolation="nearest")
                plot.add_item(array, z=j)
        else:
            for j, i in enumerate(self.surfaces):
                ofs = i.parameters.offset
                array = make.image(
                    i.array,
                    title=i.name,
                    xdata=[
                        i.parameters.job_x_origin_mm + ofs,
                        i.parameters.job_x_origin_mm + ofs +
                        i.array.shape[1] *
                        i.parameters.stage_step],
                    ydata=[
                        i.parameters.job_y_origin_mm,
                        i.parameters.job_y_origin_mm +
                        i.array.shape[0] *
                        i.parameters.stage_step],
                    interpolation="nearest")
                plot.add_item(array, z=j)
        win.exec_()
        return win

    def repeatability_test_surfaces(self,vmin=0,vmax=8):
        '''
        Method of surface dataset class
        It tests the repeatability of a set of measurements taken in the same
        area.

        Returns
        -------
        Array of the standard deviations for every point. Roughness
        repeatability
        '''

        surffla = np.array([list(i.array.flatten()) for i in self.surfaces])
        RoughnessRepetibility = np.nanstd([np.nanstd(i) for i in surffla])*1000
        stdevarr = np.nanstd(surffla,axis=0)*1000
        averagestd = round(np.nanmean(stdevarr),3)
        print('Average Measurement Repeatability: %s microns' %(averagestd))
        print('Roughness Repeatability: %s microns'%(RoughnessRepetibility))
        arr =  stdevarr.reshape(self.surfaces[0].array.shape)

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        cmap = jet
        cmap.set_bad('w', 1.)
        meshj = ax2.pcolormesh(arr, cmap = cmap)
        cm2 = plt.colorbar(meshj)
        meshj.set_clim(vmin,vmax)
        cm2.set_label('Repeatability (microns)')
        plt.show()
        return arr,RoughnessRepetibility




    def mat_export(self, mdict = None, savemask = False):
        '''
        Export in .mat format.
        '''
        import scipy.io

        self._mk_results_folder()
        if mdict == None:
            for a in self.surfaces:
                a.mat_export(savemask=savemask)
        else:
            scipy.io.savemat(path.join(
                'Results',
                self.name),
                mdict=mdict)
        self.save_log()

    def exp_csvReport(self):
        import csv
        with open('DatasetReport.csv', w) as f:
            writer = csv.writer(f)
            writer.writerow(['Sample Name',
                             'Lens name',
                             'Stage step (mm)',
                             'Calculated spot diameter (um) ',
                             'Lens repeatability (um)',
                             'CCD Freq (Hz)',
                             'X ax. speed',
                             'Y ax. speed',
                             'Laser power',
                             'Duration',
                             'Number of missing values'])
            for a in self.surfaces:
                fields = [
                    a.name,
                    a.lens.Name,
                    a.parameters.stage_step,
                    a.log['Calculated spot diameter'],
                    a.lens.Repeatability,
                    a.parameters.CCDfreq,
                    a.parameters.X_axis_vel,
                    round(a.parameters.Y_axis_vel,2),
                    a.parameters.laserpower,
                    a.parameters.time['duration'],
                    a.log['Number of missing value'],
                ]
                writer.writerow(fields)
            h,m,s = self.time_evaluation()
            writer.writerow(['',
                             '',
                             '',
                             '',
                             '',
                             '',
                             '',
                             '',
                             'Tot. Time:',
                             '%d h. %d min. %d sec.' %(h,m,s),
                             ''])


    def ROI_size_variation_effect(self,step=20,plot='show',
                                  save=False, kind='margin', unit='pt' ):

        distancefromcenter = []
        Sq=[]
        Rskw=[]
        Rkurt=[]
        surfname = []
        num_surf = len(self.surfaces)

        for n,surf in enumerate(self.surfaces):
            print(('Working on surface %s of %s:' %(n,num_surf)))
            if kind == 'margin':
                distancefromcenter1,Sq1,Rskw1, Rkurt1 = surf.margin_ROI(step,plot=False)
            else:
                distancefromcenter1,Sq1,Rskw1, Rkurt1 = surf.concentric_ROI(step,plot=False)
            if unit == 'pt':
                distancefromcenter.append(np.array(distancefromcenter1))
            elif unit == 'mm':
                distancefromcenter.append(np.array(distancefromcenter1)*surf.parameters.stage_step)
            Sq.append(Sq1)
            Rskw.append(Rskw1)
            Rkurt.append(Rkurt1)
            surfname.append(surf.name)

        if plot !=False:
            if kind == 'margin':
                xlabel = 'Margin (%s)' %(unit)
            else:
                xlabel = 'Distance from the center (%s)' %(unit)


            figx = plt.figure()
            axx = figx.add_subplot(111)
            axx.set_xlabel(xlabel)
            axx.set_ylabel('Roughness')
            axx.set_title('Roughness computed excluding different margins')


            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.set_xlabel(xlabel)
            ax2.set_ylabel('Skweness')
            ax2.set_title('Skweness computed excluding different margins')

            figt = plt.figure()
            axt = figt.add_subplot(111)
            axt.set_xlabel(xlabel)
            axt.set_ylabel('Kurtosis')
            axt.set_title('Kurtosis computed excluding different margins')

            for index,i in enumerate(self.surfaces):
                axx.plot(distancefromcenter[index],Sq[index],label = i.name)
                ax2.plot(distancefromcenter[index],Rskw[index],label = i.name)
                axt.plot(distancefromcenter[index],Rkurt[index],label =  i.name)

            #To make sure we don't have twice the same color
            colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired
            colorsx = [colormap(i) for i in np.linspace(0, 0.9, len(axx.lines))]
            colors2 = [colormap(i) for i in np.linspace(0, 0.9,len(ax2.lines))]
            colorst = [colormap(i) for i in np.linspace(0, 0.9,len(axt.lines))]
#            #linestyles = ['_', '-', '--', ':']
            #import random
            for t,j1 in enumerate(axx.lines):
                j1.set_color(colorsx[t])
               # j1.set_linestyle(random.choice[linestyles])
            for k,j2 in enumerate(ax2.lines):
                j2.set_color(colors2[k])
                #j2.set_linestyle(random.choice[linestyles])
            for z,j3 in enumerate(axt.lines):
                j3.set_color(colorst[z])
                #j3.set_linestyle(random.choice[linestyles])

            axt.legend()
            ax2.legend()
            axx.legend()
            if save:
                len_array = [len(i) for i in distancefromcenter]
                max_indx = len_array.index(max(len_array))
                self._save_data2csv(name = "Sq_vs_%s" %(kind),
                               header = list(distancefromcenter[max_indx]),
                               data = Sq,
                               surfname = surfname,
                               firstcell = kind)
                self._save_data2csv(name = "Rskw_vs_%s" %(kind),
                               header = list(distancefromcenter[max_indx]),
                               data = Rskw,
                               surfname = surfname,
                               firstcell = kind)
                self._save_data2csv(name = "Rkurt_vs_%s" %(kind),
                               header = list(distancefromcenter[max_indx]),
                               data = Rkurt,
                               surfname = surfname,
                               firstcell = kind)



        return distancefromcenter,Sq,Rskw, Rkurt


    def corners_detector(self, data = 'mask',w = 100,h = 100,
                     Kernel_1 = (15,15),
                     Kernel_2 = (200,200),
                     pad_width = 100,
                     subdivide = True,
                     method = "cv2.TM_CCORR_NORMED",
                     plot = False
                     ):
        for i in self.surfaces:
            i.corners_detector( data = data,
                                w = w,h = h,
                                Kernel_1 = Kernel_1,Kernel_2 = Kernel_2,
                                pad_width = pad_width,
                                subdivide = subdivide,
                                method = method,
                                plot = plot)


    def split_surfaces(self,padding=(0,0),dimension=100,number=3):
        newsurfdat = []
        for k in range(len(self.surfaces)):
            #we use pop to save memory
            surf = self.surfaces.pop()
            for i in range(number):
                y1 = padding[0]
                y2 = padding[0] + dimension
                x1 = padding[1] + dimension*i
                x2 = padding[1] + dimension*i + dimension
                newsurfdat.append(surf.crop(y1,y2,x1,x2,copy=True,suffix="_%s"%(i)))
        self.surfaces = newsurfdat
        self.log.append("Created %s ROI of side %s unit" %(number,dimension))


    def gwy_export(self,unit = "mm"):
        from gwyfile.objects import GwyContainer, GwyDataField, GwySIUnit
        obj = GwyContainer()

        for idx,surf in enumerate(self.surfaces):
            obj["/%s/data/title" %(idx)] = surf.name

            if  hasattr(surf.array,"mask"):
                if unit == 'mm':
                    data = GwyDataField(surf.array.data.astype(np.float))
                    data.si_unit_z = GwySIUnit(unitstr = "mm")
                if unit == 'micron' or unit == 'um':
                    data = GwyDataField(surf.array.data.astype(np.float)*1000)
                    data.si_unit_z = GwySIUnit(unitstr = "um")
                data.si_unit_xy = GwySIUnit(unitstr = "mm")
                data.xres = surf.array.shape[1]
                data.yres = surf.array.shape[0]
                data.xreal = surf.parameters.rangeX
                data.yreal = surf.parameters.rangeY

                obj["/%s/mask"%(idx)] = GwyDataField(surf.array.mask.astype(np.float))
                obj["/%s/data"%(idx)] = data
            else:
                if unit == 'mm':
                    data = GwyDataField(surf.array.astype(np.float))
                    data.si_unit_z = GwySIUnit(unitstr = "mm")
                if unit == 'micron' or unit == 'um':
                    data = GwyDataField(surf.array.astype(np.float)*1000)
                    data.si_unit_z = GwySIUnit(unitstr = "um")
                data.xres = surf.array.shape[1]
                data.yres = surf.array.shape[0]
                data.xreal = surf.parameters.rangeX
                data.yreal = surf.parameters.rangeY
                data.si_unit_xy = GwySIUnit(unitstr = "mm")
                obj["/%s/data"%(idx)] = GwyDataField(data)
            metadata = GwyContainer()
            metadata['OFFSET_mm'] = str(surf.parameters.offset)
            metadata['RANGE_X_mm'] = str(surf.parameters.rangeX)
            metadata['RANGE_Y_mm'] = str(surf.parameters.rangeY)
            metadata['PROBE_FREQUENCY_hz'] = str(surf.parameters.CCDfreq)
            metadata['JOB_X_ORIGIN_mm'] = str(surf.parameters.job_x_origin_mm)
            metadata['JOB_Y_ORIGIN_mm'] = str(surf.parameters.job_y_origin_mm)
            metadata['xAXIS_SPEED_mm/s'] = str(surf.parameters.X_axis_vel)
            metadata['STAGE_STEP_mm'] = str(surf.parameters.stage_step)
            metadata['RANGE_X_pt'] = str(surf.parameters.numrows)
            metadata['RANGE_X_pt'] = str(surf.parameters.numcols)
            metadata['LENS_NAME'] = surf.lens.Name
            metadata['LENS_PN'] = surf.lens.PN
            metadata['LASER_POWER'] = str(surf.parameters.laserpower)
            metadata['LENS_SN'] = str(surf.lens.SN)
            metadata['LENS_MINd'] = str(surf.lens.LENS_MINd)
            metadata['LENS_MAXd'] = str(surf.lens.LENS_MAXd)
            metadata['REPRODUCIBILITY'] = str(surf.lens.Reproducibility)
            metadata['REPEATABILITY'] = str(surf.lens.Repeatability)
            metadata['X_LASER_SPOT_SIZE'] = str(surf.lens.X_laser_Spot_Size)
            obj["/%s/meta"%(idx)] = metadata
    #        log = GwyContainer()
    #        log[]
        obj.tofile("%s.gwy" %(self.name))

    def hdf5_export(self, kind='minimal', compression=None,
                    toberotated =[],rotation=None):
        import h5py
        '''
        Flip to -1 if your sample is upside down.

        Rotation radians (1,2,3)
        '''

        if toberotated !=[]:
            for i in self.surfaces:
                if i.name in toberotated:
                    i.corners["Rotation"] = rotation
                    print("Flipped: %s" %(i.name))

        # create the HDF5 file
        f = h5py.File('%s_%s.h5' % (self.name, kind), "w")

        # set the attributes that are specific for the file
        f.attrs['HDF5_Version'] = h5py.version.hdf5_version
        f.attrs['h5py_version'] = h5py.version.version
        diz = self.dict_surfaces_by_name() #sort the surface by name
        for i in sorted(diz, key= diz.get):
            print(i.name)
            corners_mat = np.array([None,None])
            if i.corners["TL"] is not None and i.corners["BR"] is not None:

                #create a croner matrix
                if i.corners["Rotation"] != None:
                    if rotation == 2:
                        print("Rotating the sample")
                        sh = np.array(i.array.shape)[::-1]
                        cor = i.corners.copy()
                        i.corners["TR"] = sh - np.array(cor['BL']
                        i.corners["TR"][i.corners["TR"]<0] = 0
                        i.corners["TL"] = sh - np.array(cor['BR'])
                        i.corners["TL"][i.corners["TL"]<0] = 0
                        i.corners["BR"] = sh - np.array(cor['TL'])
                        i.corners["BR"][i.corners["BR"]<0] = 0
                        i.corners["BL"] = sh - np.array(cor['TR'])
                        i.corners["BL"][i.corners["BL"]<0] = 0
                corners_mat = np.array([[i.corners["TL"],i.corners["TR"]],
                                       [i.corners["BL"],i.corners["BR"]]])



            if kind == 'minimal':
                # Create new dataset for every surface
                datset = f.create_dataset(
                    i.name, data=i.array, compression=compression)
                datset.attrs['OFFSET_mm'] = float(i.parameters.offset)
                datset.attrs['RANGE_X_mm'] = float(i.parameters.rangeX)
                datset.attrs['RANGE_Y_mm'] = float(i.parameters.rangeY)
                datset.attrs['PROBE_FREQUENCY_hz'] = int(i.parameters.CCDfreq)
                datset.attrs['JOB_X_ORIGIN_mm'] = float(
                    i.parameters.job_x_origin_mm)
                datset.attrs['JOB_Y_ORIGIN_mm'] = float(
                    i.parameters.job_y_origin_mm)
                datset.attrs[
                    'xAXIS_SPEED_mm/s'] = float(i.parameters.X_axis_vel)
                datset.attrs['STAGE_STEP_mm'] = float(i.parameters.stage_step)
                datset.attrs['RANGE_X_pt'] = int(i.parameters.numrows)
                datset.attrs['RANGE_X_pt'] = int(i.parameters.numcols)
                datset.attrs['LENS_NAME'] = str(i.lens.Name)
                datset.attrs['LENS_PN'] = str(i.lens.PN)
                datset.attrs['LASER_POWER'] = float(i.parameters.laserpower)
                datset.attrs['LENS_SN'] = str(i.lens.SN)
                datset.attrs['LENS_MINd'] = float(i.lens.LENS_MINd)
                datset.attrs['LENS_MAXd'] = float(i.lens.LENS_MAXd)
                datset.attrs['REPRODUCIBILITY'] = float(i.lens.Reproducibility)
                datset.attrs['REPEATABILITY'] = float(i.lens.Repeatability)
                datset.attrs['X_LASER_SPOT_SIZE'] = float(
                    i.lens.X_laser_Spot_Size)
                if corners_mat.all() !=  None:
                    yi = np.min([corners_mat[0][0][1],
                                 corners_mat[0][1][1]])
                    yf = np.max([corners_mat[1][0][1],
                                 corners_mat[1][1][1]])
                    xi = np.min([corners_mat[0][0][0],
                                 corners_mat[1][0][0]])
                    xf = np.max([corners_mat[0][1][0],
                                 corners_mat[1][1][0]])
                    dataset.attrs['ROI'] = datset.regionref[yi:yf,xi:xf]

            else:
                # create a new group for every surface
                group = f.create_group(i.name)
                group.attrs['OFFSET_mm'] = float(i.parameters.offset)
                group.attrs['RANGE_X_mm'] = float(i.parameters.rangeX)
                group.attrs['RANGE_Y_mm'] = float(i.parameters.rangeY)
                group.attrs['PROBE_FREQUENCY_hz'] = int(i.parameters.CCDfreq)
                group.attrs['JOB_X_ORIGIN_mm'] = float(
                    i.parameters.job_x_origin_mm)
                group.attrs['JOB_Y_ORIGIN_mm'] = float(
                    i.parameters.job_y_origin_mm)
                group.attrs[
                    'xAXIS_SPEED_mm/s'] = float(i.parameters.X_axis_vel)
                group.attrs[
                    'yAXIS_SPEED_mm/s'] = float(i.parameters.Y_axis_vel)
                group.attrs['LASER_POWER'] = int(i.parameters.laserpower)
                group.attrs['STAGE_STEP_mm'] = float(i.parameters.stage_step)
                group.attrs['RANGE_X_pt'] = int(i.parameters.numrows)
                group.attrs['RANGE_X_pt'] = int(i.parameters.numcols)
                group.attrs['LENS_NAME'] = str(i.lens.Name)
                group.attrs['LENS_PN'] = str(i.lens.PN)
                group.attrs['LENS_SN'] = str(i.lens.SN)
                group.attrs['LENS_MINd_mm'] = str(i.lens.LENS_MINd)
                group.attrs['LENS_MAXd_mm'] = str(i.lens.LENS_MAXd)
                group.attrs['REPRODUCIBILITY'] = float(i.lens.Reproducibility)
                group.attrs['REPEATABILITY'] = float(i.lens.Repeatability)
                group.attrs['X_LASER_SPOT_SIZE_um'] = float(
                    i.lens.X_laser_Spot_Size)
                group.attrs['ACQUISITION_DATE'] = str(i.parameters.date)
                group.attrs['PROBE_FINEPOWER'] = int(i.parameters.finepower)
                group.attrs['PROBE_COURSEPOWER'] = int(i.parameters.coursepower)

                if kind == 'full' or kind == 'fullsample':
                    # if kind == 'full' all the processed data are saved
                    data = i.array
                    if i.corners['Rotation'] != None:
                        data = np.rot90(data.copy(),i.corners['Rotation'])
                    DATA = group.create_dataset(
                        '%s data' %
                        (i.name),
                        data=data,
                        compression=compression)
                    DATA.attrs['supplementary rows'] = int(
                        i.log['supplementary rows'])
                    DATA.attrs['supplementary columns'] = int(
                        i.log['supplementary columns'])
                    DATA.attrs['detected_corners'] = corners_mat
                    snr_data = i.SNR(
                            matrix=True,
                            plot=False)[1]
                    if i.corners['Rotation'] != None:
                        snr_data = np.rot90(snr_data,i.corners['Rotation'])
                    SNR = group.create_dataset(
                        '%s SNR' %
                        (i.name),
                        data= snr_data,
                        compression=compression)
                    SNR.attrs['Average SNR'] = int(i.log['SNR'])
                    if corners_mat.all() != None:
                        yi = np.min([corners_mat[0][0][1],
                                     corners_mat[0][1][1]])
                        yf = np.max([corners_mat[1][0][1],
                                     corners_mat[1][1][1]])
                        xi = np.min([corners_mat[0][0][0],
                                     corners_mat[1][0][0]])
                        xf = np.max([corners_mat[0][1][0],
                                     corners_mat[1][1][0]])
                        print(yi,yf,xi,xf)
                        DATA.attrs['ROI'] = DATA.regionref[yi:yf,xi:xf]
                        SNR.attrs['ROI'] = SNR.regionref[yi:yf,xi:xf]

                    if i.missingvaluesarray is not None:
                        missing = i.missingvaluesarray
                        if i.corners['Rotation'] != None:
                            missing = np.rot90(missing.copy(),
                                               i.corners['Rotation'])
                        MIS = group.create_dataset(
                            '%s MISSING' %
                            (i.name),
                            data=missing,
                            compression=compression)
                        MIS.attrs['NUMBER_OF_MISSINGS'] = i.log[
                            'Number of missing value']
                        MIS.attrs['MISSING_VALUE_CORRECTION'] = i.log[
                            'Missing Value correction']
                        if corners_mat.all() != None:
                            MIS.attrs['ROI'] = MIS.regionref[yi:yf,xi:xf]
                    if hasattr(i.array, 'mask'):
                        if i.array.mask.size != 1:
                            mask_data = i.array.mask
                            if i.corners['Rotation'] != None:
                                mask_data = np.rot90(mask_data.copy(),
                                                     i.corners['Rotation'])
                            MASKH = group.create_dataset(
                                '%s MASK' %
                                (i.name),
                                data=mask_data,
                                compression=compression)
                            if corners_mat.all() != None:
                                MASKH.attrs['ROI'] = MASKH.regionref[yi:yf,xi:xf]
                    if kind == 'fullsample' and i.sample_infos is not None:
                        PHOTO = group.create_dataset(
                            '%s photo' %
                            (i.name),
                            data=i.sample_infos.get_image(),
                            compression=compression)
                        PHOTO.attrs['Material'] = i.sample_infos.materials
                        PHOTO.attrs['Dimensions'] = ('%s x %s x %s cm') % (
                            i.sample_infos.width,
                            i.sample_infos.height,
                            i.sample_infos.depth)
                        PHOTO.attrs['Description'] = i.sample_infos.description

                if kind == 'raw':
                    # if kind == 'raw' the raw data are saved
                    group.create_dataset(
                        '%s data' %
                        (i.name),
                        data=np.fromfile(
                            i.path +
                            '.dist',
                            dtype=np.float32),
                        compression=compression)
                    group.create_dataset(
                        '%s SNR' %
                        (i.name),
                        data=np.fromfile(
                            i.path +
                            '.snr',
                            dtype=np.uint16),
                        compression=compression)
                    group.create_dataset(
                        '%s TAG' %
                        (i.name),
                        data=np.fromfile(
                            i.path +
                            '.tag',
                            dtype='uint8'),
                        compression=compression)


        f.close()  # Close the file

    def reference_to_corners(self,corner = 'TL',pad = 0, trim = False):
        for i in self.surfaces:
            #cut the samples so that they start at the corner
            i.crop(i.corners[corner][1]-pad,i.array.shape[0],
                   i.corners[corner][0]-pad,i.array.shape[1])
        #reshape to the same dimension
        if trim:
            shapes = np.array([i.array.shape for i in self.surfaces])
            minshape = np.min(shapes,axis = 0 )
            for i in self.surfaces:
                i.crop(0,minshape[0],0,minshape[1])

    def multi_seg_min_max(self,step,save=False,average=False,
                          xmaxlim=None,ymaxlim=None,xsquared=False,
                          export = None):
        data = {}
        count=1
        for i in self.surfaces:
            print('--------------')
            print("%s: %s of %s" %(i.name,count,len(self.surfaces)))
            print('--------------')
            data[i.name] = i.multi_segmentation(step=step,kind=6,plot=False)
            count+=1

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        fig2l = plt.figure()
        ax2l = fig2l.add_subplot(111)
        fig3l = plt.figure()
        ax3l = fig3l.add_subplot(111)
        fig4l = plt.figure()
        ax4l = fig4l.add_subplot(111)
        avymin = []
        avymax = []
        if export == 'mat':
           matfile = {}
        if export == 'csv':
            indicescsv = []
            ymincsv = []
            ymaxcsv = []
            maxrcsv = []


        for k in data:
            print('--------------')
            print(k)
            print('--------------')
            indices = np.array(sorted(data[k].keys()))

            ymin = [np.nanmean(data[k][i][0].flatten()*1000) for i in indices]
            ymax = [np.nanmean(data[k][i][1].flatten()*1000) for i in indices]
            #add also the max and min value of the whole array
            #ymin.append(np.nanmin(data[k].values()[0][0].flatten()*1000))
            #ymax.append(np.nanmax(data[k].values()[0][1].flatten()*1000))
            #indices = np.append(indices, data[k].values()[0][0].shape[0]*indices[0])
            if xsquared:
                indices = indices**2
            ax2.plot(indices,ymin,color='b',label=k)
            ax3.plot(indices,ymax,color='r',label=k)
            ax4.plot(indices,np.array(ymax)+np.array(ymin),color='r',label=k)
            ax2l.plot(indices,ymin,color='b',label=k)
            ax3l.plot(indices,ymax,color='r',label=k)
            ax4l.plot(indices,np.array(ymax)+np.array(ymin),color='r',label=k)
            if average:
                avymin.append(np.array(ymin))
                avymax.append(np.array(ymax))
            if export == 'mat':
                matfile[k+str('_indices')] = indices
                matfile[k+str('_ymin')] = ymin
                matfile[k+str('_ymax')] = ymax
                matfile[k+str('_maxr')] = np.array(ymax)+np.array(ymin)

            if export == 'csv':
                indicescsv.append(['Name']+list(indices))
                ymincsv.append([k]+ymin)
                ymaxcsv.append([k]+ymax)
                maxrcsv.append([k]+list(np.array(ymax)+np.array(ymin)))



        if average:
            avminf = np.mean(np.array(avymin),axis=0)
            avmaxf = np.mean(np.array(avymax),axis=0)
            with open('Averageymin.csv', w) as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indices)
                writer.writerow(avminf)
                writer.writerow(avmaxf)


        ax4.set_xlabel('Evaluation lenght (micron)')
        ax4.set_ylabel('Maximal  range roughness amplitude  (micron)')
        ax3.set_xlabel('Evaluation lenght (micron)')
        ax3.set_ylabel('Maximal roughness  amplitude  (micron)')
        ax2.set_xlabel('Evaluation lenght (micron)')
        ax2.set_ylabel('-Minimal roughness  amplitude  (micron)')
        ax4l.set_xlabel('Evaluation lenght (micron)')
        ax4l.set_ylabel('Maximal  range roughness amplitude  (micron)')
        ax3l.set_xlabel('Evaluation lenght (micron)')
        ax3l.set_ylabel('Maximal roughness  amplitude  (micron)')
        ax2l.set_xlabel('Evaluation lenght (micron)')
        ax2l.set_ylabel('-Minimal roughness  amplitude  (micron)')
        ax4l.set_yscale('log')
        ax3l.set_yscale('log')
        ax2l.set_yscale('log')
        ax4l.set_xscale('log')
        ax3l.set_xscale('log')
        ax2l.set_xscale('log')
        if xmaxlim != None:
            ax4.set_xlim(xmaxlim[0],xmaxlim[1])
            ax3.set_xlim(xmaxlim[0],xmaxlim[1])
            ax2.set_xlim(xmaxlim[0],xmaxlim[1])
            ax4l.set_xlim(xmaxlim[0],xmaxlim[1])
            ax3l.set_xlim(xmaxlim[0],xmaxlim[1])
            ax2l.set_xlim(xmaxlim[0],xmaxlim[1])

        if ymaxlim != None:
            ax4.set_ylim(ymaxlim[0],ymaxlim[1]*2)
            ax3.set_ylim(ymaxlim[0],ymaxlim[1])
            ax2.set_ylim(ymaxlim[0],ymaxlim[1])
            ax4l.set_ylim(ymaxlim[0],ymaxlim[1]*2)
            ax3l.set_ylim(ymaxlim[0],ymaxlim[1])
            ax2l.set_ylim(ymaxlim[0],ymaxlim[1])


        colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired
        colors = [colormap(i) for i in np.linspace(0, 0.9,len(ax2.lines))]
        for i,j in enumerate(ax2.lines):
            j.set_color(colors[i])
        for i,j in enumerate(ax3.lines):
            j.set_color(colors[i])
        for i,j in enumerate(ax4.lines):
            j.set_color(colors[i])
        for i,j in enumerate(ax2l.lines):
            j.set_color(colors[i])
        for i,j in enumerate(ax3l.lines):
            j.set_color(colors[i])
        for i,j in enumerate(ax4l.lines):
            j.set_color(colors[i])

        handles, labels = ax2.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        ax2.legend(handles, labels,prop={'size':8},loc=4)

        handles, labels = ax3.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        ax3.legend(handles, labels, prop={'size':8},loc=4)

        handles, labels = ax4.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        ax4.legend(handles, labels, prop={'size':8},loc=4)

        handles, labels = ax2l.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        ax2l.legend(handles, labels,prop={'size':8},loc=4)

        handles, labels = ax3l.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        ax3l.legend(handles, labels, prop={'size':8},loc=4)

        handles, labels = ax4l.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        ax4l.legend(handles, labels, prop={'size':8},loc=4)
        plt.tight_layout()

        self.log.append('multi_seg_min_max. Step: %s \n' %(step))
        if save:
            self._save_plot(fig2,'MinusMinimalR')
            self._save_plot(fig3,'MaxR')
            self._save_plot(fig4,'MaxRangeR')
            self._save_plot(fig2l,'MinusMinimalRLOG')
            self._save_plot(fig3l,'MaxRLOG')
            self._save_plot(fig4l,'MaxRangeRLOG')

        if export == 'mat':
           self.mat_export(matfile)
        if export == 'csv':
            with open('ymin.csv', w) as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indicescsv[0])
                for i in ymincsv:
                    writer.writerow(i)

            with open('ymax.csv', w) as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indicescsv[0])
                for i in ymaxcsv:
                    writer.writerow(i)

            with open('maxrange.csv', w) as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indicescsv[0])
                for i in maxrcsv:
                    writer.writerow(i)



    def multi_seg_stddev(self,step,save=False,average=False,
                         xlim=None,ylim=None,norm=False,export =False):
        data = {}
        count=1
        for i in self.surfaces:
            print("%s: %s of %s" %(i.name,count,len(self.surfaces)))
            data[i.name] = i.multi_segmentation(step = step,
                                                kind = 1,
                                                plot = False,
                                                norm = norm)
            count+=1

        #Normal plot
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        #Log plot
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        avstd = []
        indices = []
        if export == 'mat':
            matfile = {}

        if export == 'csv':
            indicescsv = []
            stdcsv = []
            errcsv = []


        for k in data:
            indices = np.array(sorted(data[k].keys()))
            ystd = [np.nanmean(data[k][i][0].flatten()*1000) for i in indices]
            yerr = [np.nanstd(data[k][i][0].flatten()*1000) for i in indices]
            #the std.dev of one element is 0
            ystd[0] = 0
            #ystd.append(np.nanstd(data[k][indices[0]][0].flatten()*1000))
            #indices = np.append(indices,data[k][indices[0]][0].shape[0])
            ax2.plot(indices,ystd,color='b',label=k)
            ax3.plot(indices,ystd,color='b',label=k)
            if export == 'mat':
                matfile[k+str('_indices')] = indices
                matfile[k+str('_ystd')] = ystd
            if average:
                avstd.append(np.array(ystd))
                indices = indices
            if export == 'csv':
                indicescsv.append(['Name']+list(indices))
                stdcsv.append([k]+ystd)
                errcsv.append([k]+yerr)




        if average:
            avmstd = list(np.mean(np.array(avstd),axis=0))
            with open('Averagestd.csv', w) as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indices)
                writer.writerow(avmstd)
            if export == 'mat':
                matfile[k+str('_average')] = avmstd
            return indices,avmstd



        ax2.set_xlabel('Evaluation lenght (micron)')
        ax2.set_ylabel('Roughness (micron)')
        ax3.set_xlabel('Evaluation lenght (micron)')
        ax3.set_ylabel('Roughness (micron)')
        if xlim != None:
            ax2.set_xlim(xlim[0],xlim[1])
            ax3.set_xlim(xlim[0],xlim[1])

        if ylim != None:
            ax2.set_ylim(ylim[0],ylim[1])
            ax3.set_ylim(ylim[0],ylim[1])



        colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired


        handles, labels = ax2.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        colors = [colormap(i) for i in np.linspace(0, 0.9,len(ax2.lines))]
        for i,j in enumerate(ax2.lines):
            j.set_color(colors[i])

        ax2.legend(handles, labels,prop={'size':8},loc=4)
        #FOR THE LOG PLOT
        #sort both labels and handles by labels
        handles, labels = ax3.get_legend_handles_labels()
        labels, handles = list(zip(*sorted(zip(labels, handles), key=lambda t: t[0])))
        colors = [colormap(i) for i in np.linspace(0, 0.9,len(ax3.lines))]
        for i,j in enumerate(ax3.lines):
            j.set_color(colors[i])
        ax3.legend(handles, labels,prop={'size':8},loc=4)
        ax3.set_yscale('log')
        ax3.set_xscale('log')
        plt.tight_layout()

        self.log.append('multi_seg_stddev. Step: %s' %(step))
        if save:
            self._save_plot(fig2,'Roughness')
            self._save_plot(fig3,'RoughnessLOG')

        if export == 'mat':
            self.mat_export(matfile)

        if export == 'csv':
            with open('multistd.csv', w) as f:
                    import csv
                    writer = csv.writer(f)
                    writer.writerow(indicescsv[0])
                    for i in stdcsv:
                        writer.writerow(i)
            with open('multistd_errors.csv', w) as f:
                    import csv
                    writer = csv.writer(f)
                    writer.writerow(indicescsv[0])
                    for i in errcsv:
                        writer.writerow(i)



        return indices,ystd

    def subdivide_in_multipleROIs(self,nrows,ncols,saveindexes=True):
        if self.ROIs_dataset is None:
               self.ROIs_dataset = surfacedataset()

        for num, surf in enumerate(self.surfaces):
            arr = surf.array
            h, w = arr.shape
            yicr = h/ncols
            xicr = w/nrows
            for col in range(ncols):
                for row in range(nrows):
                    y1 = yicr*row
                    y2 = yicr*row+yicr
                    x1 = xicr*col
                    x2 = xicr*col+xicr
                    suffix = "%s_%s%s" %(num,col,row)
                    roi = surf.crop(y1,y2,x1,x2,copy=True,suffix=suffix)
                    self.ROIs_dataset.surfaces.append(roi)
                    print("appended" + roi.name +suffix)
                    if saveindexes:
                        surf.ROIs_indices[suffix] = (y1,y2,x1,x2)



    def crop(self,y1,y2,x1,x2,toROIdataset=False,suffix = None,
             saveindexes = False):
        if toROIdataset:
            if self.ROIs_dataset is None:
               self.ROIs_dataset = surfacedataset()
            for i in  self.surfaces:
                surf = i.crop(y1,y2,x1,x2,copy=True,suffix=suffix)
                self.ROIs_dataset.surfaces.append(surf)
                print("appemded" + surf.name +suffix)
                if saveindexes:
                    i.ROIs_indices[suffix] = (y1,y2,x1,x2)
        else:
            for i in  self.surfaces:
                i.crop(y1,y2,x1,x2,suffix=suffix)
                if saveindexes:
                    i.ROIs_indices[suffix] = (y1,y2,x1,x2)

    def ROIinteractive_select(self,surface_index = 0,cropall= False,
                              savetoROIsDS=True,copyindextoothers=True,
                              deleteindexes=False):

        self.surfaces[surface_index].ROIinteractive_select(save='indices')
        time.sleep(1)
        for ROIkey in list(self.surfaces[surface_index].ROIs_indices.keys()):
            y1,y2,x1,x2 = self.surfaces[surface_index].ROIs_indices[ROIkey]

            if savetoROIsDS:
               self.crop(y1,y2,x1,x2,toROIdataset=True,suffix=ROIkey)

            if cropall:
                self.crop(y1,y2,x1,x2,suffix=ROIkey)


            print("done for %s" %(ROIkey))
        if copyindextoothers:
                indexes = self.surfaces[surface_index].ROIs_indices
                for i in self.surfaces:
                    i.ROIs_indices = indexes

        if deleteindexes:
            for surf in self.surfaces:
                surf.ROIs_indices = {}


    def export_spatio_temporal_profile(self,row = None,
                                       column = None,
                                       start = 0,
                                       stop = None,
                                       name = 'ouput.csv',
                                       cut_name = [0,None]):
        """

        Visualizing the results with pandas and seaborn
        data = pd.read_csv('ouput.csv',na_values = '--')
        sns.tsplot(data = data, time = 'X', value = 'Height',
        unit = 'Index',condition = 'Name')
        """
        with open(name, w) as f:
            import csv
            writer = csv.writer(f)
            writer.writerow(['X','Height','Index','Name','Date',])
            for surf_Id, surf in enumerate(self.surfaces):
                if column == None:
                    for index, val in enumerate(surf.array[row,start:stop]):
                            csv_row = [ index * surf.parameters.stage_step,
                                        val,
                                        surf_Id,
                                        surf.name[cut_name[0]:cut_name[1]],
                                        surf.parameters.date]
                            writer.writerow(csv_row)
                if row == None:
                    for index, val in enumerate(surf.array[start:stop,column]):
                            csv_row = [ index * surf.parameters.stage_step,
                                            val,
                                            surf_Id,
                                            surf.name[cut_name[0]:cut_name[1]],
                                            surf.parameters.date]
                            writer.writerow(csv_row)


    def Segmentation_Distribution(self,side,kind=[1],bins='opti',
                                  includeOL=False,normed=True,margin=None):

        figH = plt.figure()
        axb =figH.add_subplot(111)
        kinds = {1:"Sq ",2:"Ssk ",3:"Sku "}
        name = "".join([kinds[i] for i in kind])
        if bins == 'opti':
            opti = [np.linspace(0,10,100),np.linspace(-3,3,200),np.linspace(-2,2,100)]
            if includeOL:
                opti = [[-10]+i+[100] for i in opti]
            bins = [ opti[i-1] for i in kind]
        from scipy.stats.mstats import skew, kurtosis
        import csv
        with open('Segmentation_distribution_%s.csv' %(name), w) as f, open('Dataset Seg %s.csv' %(name), w) as h:
            writer = csv.writer(f)
            datasetwriter = csv.writer(h)
            writer.writerow(['Sample Name',
                             'Std. Dev.',
                             'Skweness',
                             'Kurtosis',
                             'Average'])


            for i in self.surfaces:
                if i.sample_infos != None:
                    label = i.sample_infos.label
                else:
                    label= i.name

                xdata = []
                for ki,binsx in zip(kind,bins):
                    data = i.Segmentation(side,kind=ki,margin=margin).flatten()
                    data = data[~np.isnan(data)]
                    if ki == 1:
                        data = data*1000
                    histo = axb.hist(data,bins=binsx,
                                     facecolor="r",
                                     histtype = 'step',
                                     label=label,
                                     normed=normed)
                    xdata.extend(histo[0])


                fields = [
                    i.name,
                    np.std(data),
                    skew(data),
                    kurtosis(data),
                    np.mean(data)]
                writer.writerow(fields)
                datasetwriter.writerow(xdata)



        axb.set_ylabel('Numbers of values')
        axb.set_xlabel('Roughness ( $\mu m$)')
        colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired
        colorst = [colormap(i) for i in np.linspace(0, 0.9,len(axb.patches))]
        self.log.append('Divided surfaces in pathes. Side %s' %(side))
        self.log.append('For every path has been calulated: %s' %(name))
        self.log.append('Saved Segmentation_distribution_%s.csv' %(name))
        self.log.append('Calulated Std.Dev. Skw. Kurt. and mean of the distribution')
        self.log.append('Saved the results to:Dataset Seg %s.csv' %(name))
        for t,j1 in enumerate(axb.patches):
                j1.set_color(colorst[t])
        axb.legend()
    #Utilities
    def PCA(self):
        shape = self.surfaces[0].array.shape
        nxpmatrix =[]
        for surf in self.surfaces:
            if surf.array.shape == shape:
                nxpmatrix.append(surf.array.flatten())
            else:
                print('Wrong SHAPE found!!')
        from sklearn.decomposition import PCA as sklearnPCA
        import matplotlib.cm as cm
        sklearn_pca = sklearnPCA(n_components=3)
        sklearn_transf = sklearn_pca.fit_transform(np.array(nxpmatrix))


        fig = plt.figure()
        ax = fig.add_subplot(111)

        c1exp_variance= round(sklearn_pca.explained_variance_ratio_[0]*100,2)
        c2exp_variance= round(sklearn_pca.explained_variance_ratio_[1]*100,2)
        ax.set_xlabel('First component (%s %%)' %(c1exp_variance))
        ax.set_ylabel('Second component (%s %%)' %(c2exp_variance))

        ax.scatter(sklearn_transf[:,0],sklearn_transf[:,1])
        ax.legend()
        ax.set_title('Score Plot (Correlation Matrix)')

        #Loadings

        loadings = sklearn_pca.components_

        # I've omitted the code to create ind; a list of the indexes of the
        # loadings ordered by distance from origin.

        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.scatter(*loadings, alpha=0.3, label="Loadings");


        ax3.set_title("Loading plot")
        ax3.set_xlabel("Loadings on PC1")
        ax3.set_ylabel("Loadings on PC2")
        ax3.grid()


        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        ax4.pcolormesh(loadings[0].reshape(shape))


        ax4.set_title("Loading plot")
        ax4.set_xlabel(" PC1  (pt)")
        ax4.set_ylabel("Arbitrary")


        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111)
        ax5.pcolormesh(loadings[1].reshape(shape))


        ax5.set_title("Loading plot")
        ax5.set_xlabel(" PC2  (pt)")
        ax5.set_ylabel("Arbitrary")
        self.log.append('Performed PCA')


    def align_surfaces(self,iteration=30,remove_mask=False,usemask=False):
        import imreg_dft as ird
        if usemask:
            for surf in self.surfaces[1:]:
                surf.array= ird.similarity(self.surfaces[0].array.mask,surf.array.mask)['timg']
        if remove_mask is True:
            for surf in self.surfaces:
                surf.array.mask = False
        if remove_mask =='to zero':
            for surf in self.surfaces:
                surf.array[surf.array.mask]=0
        for surf in self.surfaces[1:]:
            surf.array= ird.similarity(self.surfaces[0].array,surf.array)['timg']
        self.log.append('Surface alligment')




    def rename_with_label(self):
        '''
        Use labels instead of filenames for legends
        '''
        for sur in self.surfaces:
            sur.name = sur.sample_infos.label

    def export_LaTeX_report(self,
                            ADF=True,
                            Total=True,
                            SNR=True,
                            Timing=False,
                            cbarorientation = 'auto',
                            usemeasurementname = True):
        #Deful plt.rc params
        diz = self.dict_surfaces_by_name() #sort the surface by name
        plt.rcParams["figure.figsize"] = (20,6)
        plt.style.use('classic')
        plt.rcParams.update({'font.size':16})
        n = "\n"
        dataset_name ='dataset name'
        if self.dataset_infos != None:
            if self.dataset_infos.name != None:
                dataset_name = self.dataset_infos.name
        header = [r"\documentclass[10pt,a4paper,twoside]{article}",
        r"\usepackage[latin1]{inputenc}",
        r"\usepackage{amsmath}",
        r"\usepackage{amsfonts}",
        r"\usepackage{amssymb}",
        r"\usepackage{graphicx}",
        r"\usepackage{float}",
        r"\usepackage{hyperref}",
        r"\usepackage{lscape}",
        r"\usepackage{makecell}",
        r"\usepackage{tabularx}",
        r"\usepackage{longtable}",
        r"\renewcommand\theadalign{cb}",
        r"\renewcommand\theadfont{\bfseries}",
        r"\renewcommand\theadgape{\Gape[4pt]}",
        r"\renewcommand\cellgape{\Gape[4pt]}",
        r"\usepackage[a4paper,bindingoffset=0.2in,%",
        r"left=1in,right=1in,top=1in,bottom=1in,%",
        r"footskip=.25in]{geometry}",
        r"",
        r"\title{Report Microprofilometry:\\",
        r"%s \\ Time X \\}" %(dataset_name),
        r"\makeindex",
        r"\begin{document}",
        r"\begin{titlepage}",
        r"\centering      ",
        r"\includegraphics[width=0.7\linewidth]{Scan4recologo.png}\\[4ex]\par\vspace{1cm}",
        r"{\scshape\LARGE Universit\'a di Verona \par}      ",
        r"\vspace{1cm}    ",
        r"{\huge\bfseries Report Microprofilometry:\\     ",
        r"Bronze dataset \\ Time 3 \\\par}",
        r"\vspace{2cm}    ",
        r"{\Large\itshape Giacomo Marchioro, Claudia Daffara PhD \par}",
        r"\vfill  ",
        r"\vfill  ",
        r"\includegraphics[width=0.7\linewidth]{LOGO.png}\\[4ex]  ",
        r"% Bottom of the page    ",
        r"{\large \today\par}     ",
        r"\end{titlepage}",
        r"\newpage\null\thispagestyle{empty}\newpage",
        r"\setcounter{page}{1}",
        r"\tableofcontents",
        r"	\tableofcontents",
        r"",
        r"",
        r"	\begin{landscape}       ",
        r"\section{Summary table}",
        r"	\begin{center}  ",
        r"		\begin{longtable}{ | l| l| l| l| l| l| l| l| l| l| l | p{5cm} |}  ",
        r"			\hline  ",
        r"			 \thead{Sample Name}&  \thead{Lens name}&       ",
        r"			 \thead{Stage \\ step \\(mm)}&  ",
        r"			 \thead{Exposure \\ lenght \\ (micron)}&       ",
        r"			 \thead{Lens\\ repeatability \\(micron)}&       ",
        r"			 \thead{CCD \\Freq \\(Hz)}&     ",
        r"			 \thead{X ax.\\ Speed\\ (mm/sec)}&      ",
        r"			 \thead{Y ax.\\ Speed \\(mm/sec)}&      ",
        r"			 \thead{Laser\\ power\\ (units)}&       ",
        r"			 \thead{Duration \\(hh:mm:ss)}& ",
        r"			 \thead{Number\\ of missing \\values}   ",
        r"			 \\ \hline      ",
        r"",
        r"			\hline  ",
        r"",
        r"			\hline  "]

        body = [r"		\end{longtable}   ",
        r"",
        r"	\end{center}    ",
        r"		\end{landscape} ",
        r"		\newpage\null\thispagestyle{empty}\newpage      ",
        r"",
        r"		\twocolumn      ",
        r"		\section{Measurements}  ",
        r"",]

        for i in sorted(diz, key= diz.get):
               self._plots_save(i,ADF = ADF, Total = Total,
                                SNR = SNR, titles= False,
                                orientation = cbarorientation)

        with open('report.tex','w') as f:
           for i in header: f.write(i+n)
           for a in sorted(diz, key= diz.get):
                fields = [
                    a.name.replace("_","\\_"),
                    a.lens.Name,
                    str(a.parameters.stage_step),
                    str(a.log['Calculated spot diameter']),
                    str(a.lens.Repeatability),
                    str(a.parameters.CCDfreq),
                    str(a.parameters.X_axis_vel),
                    str(round(a.parameters.Y_axis_vel,2)),
                    str(a.parameters.laserpower),
                    str(a.parameters.time['duration']),
                    str(a.log['Number of missing value']),
                ]
                f.write(' & '.join(fields)+n)
                f.write(r"\\"+n)
                f.write(r"\hline"+n)
           h,m,s = self.time_evaluation()
           last_row = ['',
                             '',
                             '',
                             '',
                             '',
                             '',
                             '',
                             '',
                             'Tot. Time:',
                             '%d h. %d min. ' %(h,m),
                             '']
           f.write(' & '.join(last_row)+n)
           f.write(r"\\"+n)
           f.write(r"\hline"+n)
           for i in body: f.write(i+n)
           for b in sorted(diz, key= diz.get):
               #Fore each sample
                if b.sample_infos is not None:
                    description = b.sample_infos.description
                    materials = b.sample_infos.materials
                    if usemeasurementname:
                        Sample_name = b.name.replace("_","\\_")
                    else:
                        Sample_name = b.sample_infos.name.replace("_","\\_")
                else:
                    description = "Not available"
                    materials = "Not available"
                    Sample_name = b.name.replace("_","\\_")
                f.write(r"\subsection{%s}" %(Sample_name)+n)
                f.write(r"\textbf{Sample type:} %s" %(description)+n)
                f.write(r"\\      "+n)
                f.write(r"\textbf{Materials:}%s" %(materials)+n)
                f.write(r"\begin{figure}[H]      "+n)
                f.write(r"\centering      "+n)
                im = 'Results/%s im.png' %(b.name)
                f.write(r"\includegraphics[width=1\linewidth]{%s}" %(im)+n)
                f.write(r"\caption{The sample analyzed.}  "+n)
                f.write(r"\label{fig:im}    "+n)
                f.write(r"\end{figure}    "+n)
                f.write(r"\subsubsection*{Results}"+n)
                f.write(r"\begin{figure}[H]      "+n)
                f.write(r"\centering      "+n)
                adf = 'Results/%s ADF.png' %(b.name)
                f.write(r"\includegraphics[width=1\linewidth]{%s}" %(adf)+n)
                caption=""" The lines correspond to the winzoring at minus
                three sigma and plus three sigma"""
                f.write(r"\caption{The ADF function. %s}" %(caption)+n)
                f.write(r"\label{fig:adf}    "+n)
                f.write(r"\end{figure}    "+n)
                SM = 'Results/%s SM.png' %(b.name)
                f.write(r"\begin{figure}[H]      "+n)
                f.write(r"\includegraphics[width=1\linewidth]{%s}" %(SM)+n)
                f.write(r"\caption{Surface distances map. }"+n)
                f.write(r"\label{fig:adf}    "+n)
                f.write(r"\end{figure}    "+n)




                f.write(r"\subsubsection*{Acquisition parameters:}"+n)
                f.write(r"\begin{center}  "+n)
                f.write(r"\begin{tabular}{ | l| l| }      "+n)
                f.write(r"\hline  "+n)
                par1 = str(b.lens.Name)
                f.write(r"\textbf{Lens type}& %s" %(par1)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par2 = str(b.lens.X_laser_Spot_Size)
                f.write(r"\textbf{Laser spot size}& %s $\mu m$" %(par2)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par3 = str(b.lens.Repeatability)
                f.write(r"\textbf{Repeatability}& %s $\mu m$" %(par3)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par4 =str(round((b.lens.LENS_MAXd - b.lens.LENS_MINd),2))
                f.write(r"\textbf{Measurement range}& %s mm" %(par4)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par5 = str(b.parameters.laserpower)
                f.write(r"\textbf{Laser power}& %s" %(par5)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par6 =str(b.parameters.CCDfreq/1000)
                f.write(r"\textbf{CCD Freq}& %s MHz" %(par6)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par7 = str(b.parameters.X_axis_vel)
                f.write(r"\textbf{X axis speed}& %s mm/s" %(par7)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par8 = str(b.parameters.stage_step * 1000)
                f.write(r"\textbf{Stage step}& %s $\mu m$"  %(par8)+n)
                f.write(r"\\ "+n)
                f.write(r"\hline  "+n)
                f.write(r""+n)
                f.write(r""+n)
                f.write(r"\end{tabular}   "+n)
                f.write(r"\end{center}    "+n)
                #f.write(r"\newpage    "+n)
                #PROCESSING
                f.write(r"\subsubsection*{Processing:}    "+n)
                f.write(r"\begin{center}  "+n)
                f.write(r""+n)
                f.write(r"\begin{tabularx}{0.49\textwidth}{ |p{3cm}| X| }      "+n)
                f.write(r"\hline  "+n)
                par9 = b.log['Missing Value correction']
                f.write(r"\textbf{Missing Value correction}& %s" %(par9)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par10 = b.log['Best plane correction']
                f.write(r"\textbf{Form removal}& %s" %(par10) +n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par11 = b.log['Bad values filter'].replace('>',r'\textgreater')
                f.write(r"\textbf{Bad values filter}& %s " %(par11)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                f.write(r""+n)
                f.write(r""+n)
                f.write(r"\end{tabularx}   "+n)
                f.write(r"\end{center}    "+n)

                #QUALITY ASSESMENT PARAMTERS
                f.write(r"\subsubsection*{Quality assessment parameters:} "+n)
                f.write(r"\begin{center}	"	+n)
                f.write(r""+n)
                f.write(r"\begin{tabular}{ | l| l| }      "+n)
                f.write(r"\hline  "+n)
                f.write(r"\textbf{SNR}& %s" %(str(b.log['SNR']))+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par13 =str( b.log['Calculated spot diameter'])
                f.write(r"\textbf{Sampling lenght}& %s $\mu m$"
                 %(par13) + n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                par14 = str(b.log['Number of missing value'])
                f.write(r"\textbf{Missing values}& %s " %(par14)+n)
                f.write(r"\\      "+n)
                f.write(r"\hline  "+n)
                f.write(r""+n)
                par15 = str(b.log['supplementary rows'])
                f.write(r"\textbf{Deleted columns}& %s" %(par15) +n)
                f.write(r"\\      "+n)
                f.write(r"\hline	"	+n)
                f.write(r"\end{tabular}   "+n)
                f.write(r"\end{center}    "+n)
                f.write(r"\subsubsection*{Quality control:}   "+n)
                f.write(r"\begin{figure}[H]      "+n)
                f.write(r"\centering      "+n)
                s = 'Results/%s SNR.png' %(b.name)
                f.write(r"\includegraphics[width=1\linewidth]{%s}" %(s)+n)
                f.write(r"\caption{Signal to noise ratio of the sample}   "+n)
                f.write(r"\label{fig:SNR} "+n)
                f.write(r"\end{figure}    "+n)
                if Total:
                    f.write(r"\begin{figure}[H]      "+n)
                    f.write(r"\centering      "+n)
                    s1 = 'Results/%s Total.png' %(b.name)
                    f.write(r"\includegraphics[width=1\linewidth]{%s}" %(s1)+n)
                    f.write(r"\caption{Total of the sample}   "+n)
                    f.write(r"\label{fig:total}    "+n)
                    f.write(r"\end{figure}    "+n)
                if Timing:
                    f.write(r"%Timing diagram "+n)
                    f.write(r"\begin{figure}[H]      "+n)
                    f.write(r"\centering      "+n)
                    s2 = 'Results/%s SM.png' %(b.name)
                    f.write(r"\includegraphics[width=1\linewidth]{%s}"%(s2)+n)
                    f.write(r"\caption{Timing diagram of the sample}  "+n)
                    f.write(r"\label{fig:TM}    "+n)
                    f.write(r"\end{figure}    "+n)
                f.write(r""+n)
                f.write(r"\newpage      "+n)
           f.write(r"\end{document}"+n)



    def _plots_save(self, surf,
                    flipyax=None,
                    ADF=True,
                    Total=True,
                    SNR=True,
                    titles = True,
                    orientation = 'auto'):
        '''
        exp_appyreport export report with data directly from the project. The
        template used is called template.odt. The variables can be acessed
        using conditional fiedls Ctrl+F2 from the .odt document and placing 'true'
        in the 'Condition' field and the name of the variable in the 'Then' field.
        All the local variabels can be accessed e.g. self.name nested variables
        must be reassign e.g. stage_step = self.parameters.stage_step.

        '''
        import gc
        print('Plotting...')
        # Save reslts

        title = True if titles else False
        plotF= surf.plot(mode=2,plot=False,title = title,
                         orientation = orientation)
        try:
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
            self._save_plot(plotF,'%s SM' %(surf.name))


        except MemoryError:
            print('Size exceed memory limits,diluting data...')
            plotF = self.plot(mode=2, plot=False, unit='micron',dilution=2,
                              title = title, orientation = orientation)
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
            self._save_plot(plotF,'%s SM' %(surf.name))

        plt.close()
        gc.collect()
        print('ADF...')
        # ADF Image
        if ADF:
            ADFv = surf.ADF(plot=False,title = title, vminx = '-3sigma',
                            vmaxx = '+3sigma')
            self._save_plot(ADFv,'%s ADF' %(surf.name))
            plt.close()



        if surf.sample_infos is not None:

            print('Sample image..')
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.imshow(surf.sample_infos.get_image())
            ax2.tick_params(which='both', bottom='off', top='off',
                            labelbottom='off', labelleft='off',
                            left='off', right='off')
            self._save_plot(plt,'%s im' %(surf.name))
            plt.close()

        # SNR image
        print('SNR...')
        if SNR:
            plotSNR = surf.SNR(plot=True, cmap = 'd',title = title,
                               orientation = orientation)[0]
            try:
                if flipyax == -1:
                    plt.gca().invert_yaxis()
                    plt.gca().invert_xaxis()
                self._save_plot(plotSNR,'%s SNR' %(surf.name))

            except MemoryError:
                print('Size exceed memory limits,diluting data...')
                plotSNR = surf.SNR(plot=True,dilution=2,
                                   title = title, orientation = orientation)[0]
                if flipyax == -1:
                    plt.gca().invert_yaxis()

                self._save_plot(plotSNR,'%s SNR' %(surf.name))
            plt.close()

            gc.collect()

        if Total:
            plotTotal = surf.Total(plot=True, cmap ='d',
                                   title = title, orientation = orientation)[0]
            try:
                if flipyax == -1:
                    plt.gca().invert_yaxis()
                    plt.gca().invert_xaxis()
                self._save_plot(plotTotal,'%s Total' %(surf.name))

            except MemoryError:
                print('Size exceed memory limits,diluting data...')
                plotTotal = surf.Total(plot=True,dilution=2,
                                       title = title, orientation = orientation)[0]
                if flipyax == -1:
                    plt.gca().invert_yaxis()

                self._save_plot(plotTotal,'%s Total' %(surf.name))
            plt.close()

            gc.collect()

    def get_QualityControl(self):
        for i in  self.surfaces:
            i._get_QualityControl()
    def _mk_results_folder(self):
         if not path.exists('Results'):
            makedirs('Results')


    def _save_plot(self,fig,name,dpi = 100):
        self._mk_results_folder()
        fig.savefig(path.join('Results',name), dpi = dpi)

    def _save_data2csv(self,name,header,data,surfname,firstcell):
        '''
         Save data using this format
         str  firstcell | [header ... ]
         [srufname ][0] | [data ][0]
         [srufname ][1] | [data ][1]
         [srufname ][2] | [data ][2]
        '''
        with open('%s.csv' %(name), w) as f:
            import csv
            writer = csv.writer(f)
            writer.writerow([firstcell]+header)
            for j in range(len(data)):
                writer.writerow([surfname[j]]+data[j])
