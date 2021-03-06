# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a script to elaborate data form ths profilometer.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import jet, Greys
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from os import listdir, getcwd, path, makedirs,remove

import time

class Lens:

    def __init__(
            self,
            Name,
            PN,
            SN,
            LENS_MINd,
            LENS_MAXd,
            Reproducibility,
            Repeatability,
            X_laser_Spot_Size,
            angular_coverage):
        self.Name = Name
        self.PN = PN
        self.SN = SN
        self.LENS_MINd = LENS_MINd
        self.LENS_MAXd = LENS_MAXd
        self.Reproducibility = Reproducibility
        self.Repeatability = Repeatability
        self.X_laser_Spot_Size = X_laser_Spot_Size
        self.angular_coverage= angular_coverage


l50 = Lens('LENS 50mm', '3Z81050', 22263, 40, 48, 1, 0.1, 37,170)
l50bis = Lens('LENS 50mm', '3Z81050', 30225, 42, 46, 1, 0.1, 37,170)
l75 = Lens('LENS 75 mm', '3Z81075', 30225, 61, 79, 2, 0.3, 47,170)
l100 = Lens('LENS 100 mm', '3Z81100', 30225, 77.5, 112.5, 4, 0.5, 63,170)
l200 = Lens('LENS 200 mm', '3Z82007', 30225, 137.5, 262.5, 25, 3, 105,170)
l50t = Lens('LENS=50mm T', '3Z81050T', 30246, 41, 43, 0.5, 0.1, 19,150)
l25a = Lens('LENS=25mm A', '3Z8103', 3024630225, 17.95, 18.55, 0.2, 0.06, 12,150)
l25N = Lens('LENS=25mm N', '3Z79030', 000000, 15.9, 16.9, 0.06, 0.06,5,5)
l50N = Lens('LENS=50mm N','3Z799050',000000,37.5,42.5,0,0,16,3)
l75N = Lens('LENS=75mm N','3Z799050',000000,60.5,69.5,0,0,25,1.5)

diclens = {
    "l50": l50,
    "l50bis": l50bis,
    "l75": l75,
    "l100": l100,
    "l200": l200,
    "l50H": l50t,
    "l25A": l25a,
    "l25N": l25N,
    "l50N": l50N,
    "l75N": l75N,
    }


class ns:
    '''
    New surface class (ns):
    This class create an object of the micorprofilometer it allows different method
    to be used with the dot notation. The first thing to do is to create a new
    surface. This class automatically inizialize the instances with a method.

    NOTE this module is design to scan over the X axis most of the feature works
    only with this set-up.

    '''

    class parameterslist:

        def __init__(self):
            self.offset = None
            self.rangeX = None
            self.rangeY = None
            self.CCDfreq = None
            self.job_x_origin_mm = None
            self.job_y_origin_mm = None
            self.X_axis_vel = None
            self.Y_axis_vel = None
            self.stage_step = None
            self.colorbar_label = None
            self.xlabel = None
            self.report = False
            self.numrows = None
            self.numcols = None
            self.date = None
            self.laserpower = None
            self.finepower = None
            self.coursepower = None
            self.time = {'start': None, 'finish': None, 'duration': None}
            self.plane_subtraction_coef = None 

    def __init__(self, path, lens=None,init=True):
        self.path = path
        self.name = path
        self.array = np.array([])
        self.missingvaluesarray = None
        self.parameters = self.parameterslist()
        self.lens = lens
        self.sample_infos = None
        self.header = {}
        self.ROI = None
        self.corners = {"TL":None,"TR":None,"BR":None,"BL":None,"Rotation":None}
        self.ROIs_indices = {}
        self.special_ID = None
        self.QCresult = None
        self.logQC = {
            'plot-quality':None,
            'workingdistance':None,
            'SNR-QC':None,
        }
        self.log = {
            'supplementary rows': '0',
            'supplementary columns': '0',
            'removed col': 'None',
            'method': '',
            'Number of missing value': 'n.d',
            'Missing Value correction': 'Not performed',
            'SNR': 'n.d',
            'Best plane correction': 'Not performed',
            'Bad values filter': 'Not performed',
            'Processing': '',
            'Calculated spot diameter': 'n.d'}
        if init:
            try:
                self.inizialize()                
            except IOError:
                print 'No files found!'


#    def __add__(self, other):
#         newsurface = self
#         newsurface.array=self.array+other.array
#         return newsurface
#
#    def __sub__(self,other):
#        newsurface = self
#        newsurface.array=self.array-other.array
#        return newsurface
    def acoustic_init(self,path=True):
        import scipy.io  as scio
        self.array = scio.loadmat(self.name)['data']
        base = self.name.split('_TOF.mat')[0].replace('-','_')
        self.array = np.ma.masked_invalid(self.array)
        if path:
            base = getcwd() + '\\' +base
    
        with open(base + '.acs_Info.txt') as f:
            for i in f:
                g = i.split(':')
                if len(g) > 2:
                    self.header[g[0]] = g[1:]
                else:
                    self.header[g[0]] = g[1]
        self.parameters.job_x_origin_mm  = 0
        self.parameters.job_y_origin_mm = 0
        self.parameters.rangeX = float(self.header['Scanning Area x'].split()[0])
        self.parameters.rangeY = float(self.header['Scanning Area y'].split()[0])
        self.parameters.numcols = float(self.header['Columns'].split()[0])
        self.parameters.numrows = float(self.header['Rows'].split()[0])
        self.parameters.stage_step = float(self.header['Scanning Step x'].split()[0])
        self.parameters.offset = 0 
        

        
        return g
        
    def inizialize(self):
        '''
        This is a inzizlizing fucntion design to retrieve data form the header
        file of the microprofilometer. It returns all the parameters useful to
        perform further analysis of the data.
        A 2D array of raw data, unflipped and unmasked is formed at the end
        of the process and stored in self.array. This will be the base for
        further processing.
        '''
        self.path = self.path.replace(
            '.', '')  # remove any superfuls full stop
        self.name = self.path.split('\\')[-1]
        with open(self.path + '.hdr') as f:
            for i in f:
                g = i.split()
                if len(g) > 2:
                    self.header[g[0]] = g[1:]
                else:
                    self.header[g[0]] = g[1]

        print "WORKING ON:"
        print self.name
        #********************************************
        #**************HEADER VERSION 1**************
        #********************************************
        if 'RANGE_X' in self.header:
            import time
            print 'Using VERSION 1 of the header'
            self.parameters.time['start'] = self.header['JOB_TIME'][1]
            self.parameters.time['finish'] = time.ctime(
                path.getmtime(self.path + '.dist')).split()[3]
            from datetime import datetime
            s1 = self.parameters.time['start']
            s2 = self.parameters.time['finish']
            FMT = '%H:%M:%S'
            td = datetime.strptime(s2, FMT) - datetime.strptime(s1, FMT)
            self.parameters.time['duration'] = '%s:%s:%s' % (
                td.seconds // 3600, (td.seconds // 60) % 60, td.seconds % 60)
            self.parameters.stage_step = float(
                self.header['PIXEL_SIZE'].replace(
                    ',', '.'))  # get stage_step and change
            # comma with fullstop
            try:
                self.parameters.offset = float(
                    self.header['OFFSET'].replace(',', '.'))
            except:
                self.parameters.offset = input(
                    'Please define offset manually: ')

            self.parameters.rangeX = float(self.header['RANGE_X'].replace(
                ',', '.')) - self.parameters.offset * 2  # retrieve Range X minus offset
            # from .hdr fil
            self.parameters.rangeY = float(
                self.header['RANGE_Y'].replace(
                    ',', '.'))  # Retrive Range Y from .hdr file
            if self.lens is None:
                try:
                    self.lens = diclens[
                        'l' + self.header['LENS_ID'].strip('\x00')]
                    self.lens.LENS_MINd = float(self.header['LENS_MINd'])
                    self.lens.LENS_MAXd = float(self.header['LENS_MAXd'])
                except IndexError:
                    lens = raw_input(
                        "Write lens name: l50,l75,l100,l200,l50H,l25A:  ")
                    self.lens = diclens[lens]

            self.parameters.job_x_origin_mm = float(self.header['JOB_ORIGIN'][0].replace(
                ',', '.'))
            self.parameters.job_y_origin_mm = float(self.header['JOB_ORIGIN'][1].replace(
                ',', '.'))
            self.parameters.CCDfreq = float(self.header['PROBE_FREQUENCY'])
            self.parameters.X_axis_vel = self.header[
                'xAXIS_SPEED'].replace(',', '.')
            self.parameters.Y_axis_vel = float(self.header[
                'yAXIS_SPEED'].replace(',', '.'))
            self.parameters.laserpower = int(
                self.header['PROBE_POWER'].replace(',', '.'))
            if int(self.header['JOB_ID']) != 1:
                from LabDBwrapper import Sample
                self.sample_infos = Sample(int(self.header['JOB_ID']),
                                           r"C:\Users\OPdaTe\Documents\LabDatabase\SamplesDB")

        firstmeasurment = 1
        #************************************
        #******VERSION 2 OF THE HEADER ******
        #************************************

        if 'RANGE_X_pt' in self.header:
            import time
            self.parameters.date = self.header['JOB_TIME'][0]
            self.parameters.time['start'] = self.header['JOB_TIME'][1]
            self.parameters.time['finish'] = time.ctime(
                path.getmtime(self.path + '.dist')).split()[3]
            from datetime import datetime
            s1 = self.parameters.time['start']
            s2 = self.parameters.time['finish']
            FMT = '%H:%M:%S'
            td = datetime.strptime(s2, FMT) - datetime.strptime(s1, FMT)
            self.parameters.time['duration'] =  '%s:%s:%s' % (
                td.seconds // 3600, (td.seconds // 60) % 60, td.seconds % 60)
            if int(self.header['JOB_ID']) != 1:
                from LabDBwrapper import Sample
                self.sample_infos = Sample(int(self.header['JOB_ID']),
                                           r"C:\Users\OPdaTe\Documents\LabDatabase\SamplesDB")
            # stage_step
            self.parameters.stage_step = float(
                self.header['JOB_RESOLUTION_mm'].replace(',', '.'))
            # laser power
            self.parameters.laserpower = int(
                self.header['PROBE_POWER'].replace(',', '.'))
            try:
                self.parameters.offset = float(
                    self.header['OFFSET_mm'].replace(',', '.'))
            except:
                self.parameters.offset = input(
                    'Please define offset manually: ')

            self.parameters.rangeX = float(self.header['RANGE_X_mm'].replace(
                ',', '.')) - self.parameters.offset * 2  # retrieve Range X minus offset
            # from .hdr fil
            if self.header['RANGE_Y_mm'] == 'Null':
                self.parameters.rangeY = float(
                    self.header['RANGE_Y_pt']) * self.parameters.stage_step
                firstmeasurment = 0
            else:
                self.parameters.rangeY = float(
                    self.header['RANGE_Y_mm'].replace(
                        ',', '.'))  # Retrive Range Y from .hdr file
            if self.lens is None:
                try:
                    self.lens = diclens[
                        'l' + self.header['LENS_ID'].strip('\x00')]
                    self.lens.LENS_MINd = float(self.header['LENS_MINd'])
                    self.lens.LENS_MAXd = float(self.header['LENS_MAXd'])
                except IndexError:
                    lens = raw_input(
                        "Write lens name: l50,l75,l100,l200,l50H,l25A:  ")
                    self.lens = diclens[lens]

            self.parameters.job_x_origin_mm = float(
                self.header['JOB_X_ORIGIN_mm'])
            self.parameters.job_y_origin_mm = float(
                self.header['JOB_Y_ORIGIN_mm'])
            self.parameters.CCDfreq = float(self.header['PROBE_FREQUENCY_hz'])
            self.parameters.X_axis_vel = float(
                self.header['xAXIS_SPEED_mm/s'].replace(',', '.'))
            self.parameters.Y_axis_vel = float(
                self.header['yAXIS_SPEED_mm/s'].replace(',', '.'))
            self.parameters.finepower = float(
                self.header['PROBE_FINEPOWER'].replace(',', '.'))
            self.parameters.coursepower = float(
                self.header['PROBE_COURSEPOWER'].replace(',', '.'))
        if 'distance_x	' in self.header:
            import time
#            self.parameters.date = self.header['JOB_TIME'][0]
#            self.parameters.time['start'] = self.header['JOB_TIME'][1]
#            self.parameters.time['finish'] = time.ctime(
#                path.getmtime(self.path + '.dist')).split()[3]
#            from datetime import datetime
#            s1 = self.parameters.time['start']
#            s2 = self.parameters.time['finish']
#            FMT = '%H:%M:%S'
#            td = datetime.strptime(s2, FMT) - datetime.strptime(s1, FMT)
#            self.parameters.time['duration'] =  '%s:%s:%s' % (
#                td.seconds // 3600, (td.seconds // 60) % 60, td.seconds % 60)
#            if int(self.header['JOB_ID']) != 1:
#                from LabDBwrapper import Sample
#                self.sample_infos = Sample(int(self.header['JOB_ID']),
#                                           r"C:\Users\OPdaTe\Documents\LabDatabase\SamplesDB")
            # stage_step
            self.parameters.stage_step = float(
                self.header['sampling_mm'])
            # laser power
            self.parameters.laserpower = int(
                self.header['course_Pwr'])
            try:
                self.parameters.offset = float(
                    self.header['OFFSET_mm'])
            except:
                self.parameters.offset = input(
                    'Please define offset manually: ')

            self.parameters.rangeX = float(self.header['distance_x']) - self.parameters.offset * 2  # retrieve Range X minus offset
            # from .hdr fil
#            if self.header['distance_y'] == 'Null':
#                self.parameters.rangeY = float(
#                    self.header['distance_y']) * self.parameters.stage_step
#                firstmeasurment = 0
#            else:
            self.parameters.rangeY = float(
                    self.header['distance_y'])  # Retrive Range Y from .hdr file
#            if self.lens is None:
#                try:
#                    self.lens = diclens[
#                        'l' + self.header['LENS_ID'].strip('\x00')]
#                    self.lens.LENS_MINd = float(self.header['LENS_MINd'])
#                    self.lens.LENS_MAXd = float(self.header['LENS_MAXd'])
#                except IndexError:
#                    lens = raw_input(
#                        "Write lens name: l50,l75,l100,l200,l50H,l25A:  ")
#                    self.lens = diclens[lens]

            self.parameters.job_x_origin_mm = float(
                self.header['origin_x'])
            self.parameters.job_y_origin_mm = float(
                self.header['origin_y'])
            self.parameters.CCDfreq = float(self.header['ccd_Freq'])
            self.parameters.X_axis_vel = float(
                self.header['speed_x_mm/s'])
            self.parameters.Y_axis_vel = float(
                self.header['speed_y_mm/s'])
            self.parameters.finepower = float(
                self.header['fine_Pwr'])
            self.parameters.coursepower = float(
                self.header['course_Pwr'])
            
        self.log = {
            'supplementary rows': '0',
            'supplementary columns': '0',
            'removed col': 'None',
            'method': '',
            'Number of missing value': 'n.d',
            'Missing Value correction': 'Not performed',
            'SNR': 'n.d',
            'Best plane correction': 'Not performed',
            'Bad values filter': 'Not performed',
            'Processing': '',
            'Calculated spot diameter': 'n.d'}
        self.parameters.numrows = int(
            self.parameters.rangeY / self.parameters.stage_step) + firstmeasurment
        # The 1 is because the first mesure is at 0
        print self.parameters.rangeX
        self.parameters.numcols = int(
            self.parameters.rangeX / self.parameters.stage_step) + 1
        print 'Num. rows: %s, Num. Cols: %s' % (self.parameters.numrows, self.parameters.numcols)

        # LOADING THE DATA
        arr = np.fromfile(self.path + '.dist', dtype=np.float32)
        print 'Size expected: %s, Size retrieved: %s' % (self.parameters.numrows * self.parameters.numcols, arr.size)
        self.parameters.colorbar_label = '($mm$) Repeatability=%s $\mu m$' % (
            self.lens.Repeatability)

        if self.parameters.numrows * self.parameters.numcols != arr.size:
            print "WARNING! SIZES DO NOT MATCH!"
            diff = self.parameters.numrows * self.parameters.numcols - arr.size
            print 'Difference: %s' % (diff)
            print 'Maybe rows: %s' % (diff / self.parameters.numrows)
            print 'Maybe cols: %s' % (diff / self.parameters.numcols)
            try:
                b = arr.reshape(-1, self.parameters.numcols)
            except ValueError:
                g = input('Add cols manually:')
                b = arr.reshape(-1, self.parameters.numcols + g)
            srx = b.shape[0] - self.parameters.numrows
            scx = b.shape[1] - self.parameters.numcols
            print " I am trying to fix it passing %s for rows and %s for coloumns" % (srx, scx)            
            self.parameters.numrows += srx
            self.parameters.numcols += scx
            self.log['supplementary rows'] = srx
            self.log['supplementary columns'] = scx
            print "Now I am correcting the x-y range"
            self.parameters.rangeY += srx*self.parameters.stage_step
            self.parameters.rangeX += scx*self.parameters.stage_step
            
        else:
            b = arr.reshape(-1, self.parameters.numcols)
        self.array = b
        
        


    def Total(self, matrix=False, plot='show', 
            cartesian=True, absolute=False, replacewith=0,dilution=False,
            QC=False,cmap='d',title=True, orientation = 'auto'
            ):
        '''
        This function retrive data from the .snr file where the
        signal to noise ratio for every measurement is stored.
        It allows to the quality of the measurement.

        ------------------
         Keyword arguments:
             matrix --compute the matrix (default False)
             plot -- return a plot of the data (default 'show') to show the data
             use the argument plot='show' or plt.show()


        -----------------
        Returns:
            plt -- plot object
            b -- 2d array where data is reshaped and flipped if matrix = True
            otherwise return None

        '''

        a = np.fromfile(self.path + '.total', dtype=np.uint16)
        b = 'None'
        figs = None
        if orientation == 'auto':
            if self.array.shape[0] < self.array.shape[1]:
                orientation = 'horizontal'
            else:
                orientation = 'vertical'

        if plot:
            matrix = True

        if matrix:

            if self.missingvaluesarray is not None:
                # See missing value function for a better expanation of how
                # the correction is perfromed
                b1 = a.reshape(-1, self.parameters.numcols)
                print 'Now correcting the missing value...'
                repariedmatrix = []
                for i, j in zip(b1, self.missingvaluesarray):
                    leng = np.count_nonzero(j)
                    indexes = np.nonzero(j)
                    todelete = indexes[0][-leng / 2:]
                    corr = np.delete(i, todelete)
                    toinsert = indexes[0][:-leng / 2]
                    for d in toinsert:
                        corr = np.insert(corr, d, replacewith)

                    if len(corr) < len(i):
                        corr = np.append(
                            corr, np.full(
                                len(i) - len(corr), replacewith))
                    if len(corr) > len(i):
                        sl = (len(i) - len(corr))
                        corr = corr[:sl]

                    repariedmatrix.append(corr)
                print 'Reforming the final matrix...'
                b = np.array(repariedmatrix)
            else:
                b = a.reshape(-1, self.parameters.numcols)

            rowtoflip = range(1, self.parameters.numrows, 2)
            b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
            if plot:
                figs = plt.figure()
                axs = figs.add_subplot(111, aspect='equal')
                if cartesian:
                    b = np.flipud(b)
                if dilution is not False:
                    b = b[::dilution,::dilution]
       
                if cmap == 'c':
                    mesh = axs.pcolormesh(b)
                    cbar = figs.colorbar(mesh, orientation=orientation)                
                
                if cmap == 'd':
                    bounds = np.array([900,1200,16000,18000])
                    import matplotlib.colors as colors
                    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
                    mesh = axs.pcolormesh(b, norm=norm, cmap='rainbow')
                    mesh.cmap.set_over('white')
                    mesh.cmap.set_under('black')
                    cbar = figs.colorbar(mesh, 
                                        ax=axs,
                                        extend='both',
                                        orientation=orientation)
                    cbar.set_ticks([900,1050,1200,8600,16000,17000,18000])
                    cbar.set_ticklabels([900,'aceptable',1200,'optimal',16000,'aceptable',18000])                
                
                
                
                
                cbar.set_label('Total (Raw value)')

                if absolute:
                    start, end = axs.get_xlim()
                    xticks = ticker.FuncFormatter(
                        lambda x,
                        pos: '{0:g}'.format(
                            x *
                            self.parameters.stage_step +
                            float(
                                self.parameters.job_x_origin_mm) +
                            self.parameters.offset))
                    yticks = ticker.FuncFormatter(
                        lambda x,
                        pos: '{0:g}'.format(
                            x *
                            self.parameters.stage_step +
                            float(
                                self.parameters.job_y_origin_mm)))
                    axs.xaxis.set_ticks(
                        np.arange(
                            start,
                            end,
                            1 /
                            self.parameters.stage_step))
                    axs.yaxis.set_ticks(
                        np.arange(
                            start,
                            end,
                            1 /
                            self.parameters.stage_step))
                    axs.xaxis.set_major_formatter(xticks)
                    axs.yaxis.set_major_formatter(yticks)
                    axs.set_xlabel(
                        '(mm) stage step=%i $\mu m$ ' %
                        (self.parameters.stage_step * 1000))
                    axs.set_ylabel('(mm)')

                axs.set_xlim(0, b.shape[1])
                axs.set_ylim(0, b.shape[0])
                if cartesian == False:
                    axs.invert_yaxis()
                if title:
                    axs.set_title('Total plot %s' % (self.name))
                
                if QC:
                    from matplotlib.widgets import Button
                    plt.subplots_adjust(bottom=0.2)
                    def update(status):   
                            self.logQC['Total-QC'] = status
                            plt.close(figs)
        
                    class Index:
                        ind = 0
                        def passed(self, event):
                            update('Passed')
                            
                    
                        def failed(self, event):
                            update('Failed')
                   
                    
                    callback = Index()
                    axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
                    axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
                    bnext = Button(axnext, 'PASSED', color='g',hovercolor='y')
                    bnext.on_clicked(callback.passed)
                    bprev = Button(axprev, 'FAILED',color='r',hovercolor='y')
                    bprev.on_clicked(callback.failed)
                if plot == 'show':
                    plt.show()

        self.log['Total'] = round(np.mean(a), 2)
        print 'Average Total: %s' % (self.log['Total'])
        return figs, b
        
    def SNR(self, matrix=False, plot='show', 
            cartesian=True, absolute=False, replacewith=0,dilution=False,
            display='percentage',QC = False,cmap = 'c',title=True,
            orientation = 'auto' ):
        '''
        This function retrive data from the .snr file where the
        signal to noise ratio for every measurement is stored.
        It allows to the quality of the measurement.

        ------------------
         Keyword arguments:
             matrix --compute the matrix (default False)
             plot -- return a plot of the data (default 'show') to show the data
             use the argument plot='show' or plt.show()


        -----------------
        Returns:
            plt -- plot object
            b -- 2d array where data is reshaped and flipped if matrix = True
            otherwise return None

        '''

        a = np.fromfile(self.path + '.snr', dtype=np.uint16)
        b = 'None'
        figs = None
        if orientation == 'auto':
            if self.array.shape[0] < self.array.shape[1]:
                orientation = 'horizontal'
            else:
                orientation = 'vertical'

        if plot:
            matrix = True

        if matrix:

            if self.missingvaluesarray is not None:
                # See missing value function for a better expanation of how
                # the correction is perfromed
                b1 = a.reshape(-1, self.parameters.numcols)
                print 'Now correcting the missing value...'
                repariedmatrix = []
                for i, j in zip(b1, self.missingvaluesarray):
                    leng = np.count_nonzero(j)
                    indexes = np.nonzero(j)
                    todelete = indexes[0][-leng / 2:]
                    corr = np.delete(i, todelete)
                    toinsert = indexes[0][:-leng / 2]
                    for d in toinsert:
                        corr = np.insert(corr, d, replacewith)

                    if len(corr) < len(i):
                        corr = np.append(
                            corr, np.full(
                                len(i) - len(corr), replacewith))
                    if len(corr) > len(i):
                        sl = (len(i) - len(corr))
                        corr = corr[:sl]

                    repariedmatrix.append(corr)
                print 'Reforming the final matrix...'
                b = np.array(repariedmatrix)
            else:
                b = a.reshape(-1, self.parameters.numcols)

            rowtoflip = range(1, self.parameters.numrows, 2)
            b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
            if plot:
                figs = plt.figure()
                axs = figs.add_subplot(111, aspect='equal')
                if cartesian:
                    b = np.flipud(b)
                if dilution is not False:
                    b = b[::dilution,::dilution]
                    
                if display == 'percentage':
                    b = b/10.24
                    bounds = np.array([30,50,100])
                    tiks = [30,40,50,75,100]
                    tiks_lab = [30,'aceptable',50,'accurate',100]
                    labelKJ = 'SNR (%)'
                else:
                    mesh = axs.pcolormesh(b)                   
                    bounds = np.array([307,512,1024])
                    labelKJ = 'SNR (Raw value)'
                    tiks = [307,409,512,768,1024]
                    tiks_lab = [307,'aceptable',512,'accurate',1024]
                
                if cmap == 'c':
                    mesh = axs.pcolormesh(b)
                    cbar = figs.colorbar(mesh, orientation=orientation)                
                
                if cmap == 'd':
                    
                    import matplotlib.colors as colors
                    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
                    mesh = axs.pcolormesh(b, norm=norm, cmap='brg')
                    mesh.cmap.set_under('black')
                    cbar = figs.colorbar(mesh, 
                                        ax=axs,
                                        extend='min',
                                        orientation=orientation)
                    cbar.set_ticks(tiks)
                    cbar.set_ticklabels(tiks_lab)                
                
                cbar.set_label(labelKJ)
                
                if absolute:
                    start, end = axs.get_xlim()
                    xticks = ticker.FuncFormatter(
                        lambda x,
                        pos: '{0:g}'.format(
                            x *
                            self.parameters.stage_step +
                            float(
                                self.parameters.job_x_origin_mm) +
                            self.parameters.offset))
                    yticks = ticker.FuncFormatter(
                        lambda x,
                        pos: '{0:g}'.format(
                            x *
                            self.parameters.stage_step +
                            float(
                                self.parameters.job_y_origin_mm)))
                    axs.xaxis.set_ticks(
                        np.arange(
                            start,
                            end,
                            1 /
                            self.parameters.stage_step))
                    axs.yaxis.set_ticks(
                        np.arange(
                            start,
                            end,
                            1 /
                            self.parameters.stage_step))
                    axs.xaxis.set_major_formatter(xticks)
                    axs.yaxis.set_major_formatter(yticks)
                    axs.set_xlabel(
                        '(mm) stage step=%i $\mu m$ ' %
                        (self.parameters.stage_step * 1000))
                    axs.set_ylabel('(mm)')

                axs.set_xlim(0, b.shape[1])
                axs.set_ylim(0, b.shape[0])
                if cartesian == False:
                    axs.invert_yaxis()
                
                if title:
                    axs.set_title('SNR plot %s' % (self.name))
                
                if QC:
                    from matplotlib.widgets import Button
                    plt.subplots_adjust(bottom=0.2)
                    def update(status):
                            self.logQC['SNR-QC'] = status
                            plt.close(figs)
        
                    class Index:
                        ind = 0
                        def passed(self, event):
                            update('Passed')
                            
                    
                        def failed(self, event):
                            update('Failed')
                   
                    
                    callback = Index()
                    axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
                    axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
                    bnext = Button(axnext, 'PASSED', color='g',hovercolor='y')
                    bnext.on_clicked(callback.passed)
                    bprev = Button(axprev, 'FAILED',color='r',hovercolor='y')
                    bprev.on_clicked(callback.failed)
                if plot == 'show':
                    plt.show()

        self.log['SNR'] = round(np.mean(a), 2)
        print 'Average SNR: %s' % (self.log['SNR'])
        return figs, b
    
    def find_best_sphere(self):
        for line in self.array:
            pass
            
        
    def diagnose(self, b, numcols, numrows, v=False):
        '''
        Method of ns class.
        This function is used from other methods to diagnose errors.
        It perfroms a series of test over the data.

        Parameters
        ----------
        b : np.array
                 matrix
        numcols : int
                 number of columns caclulated
        numrows : int
                 number of rows calculated
        v : Boole
                 activate verbose mode (default is False)

        Returns
        -------
        Prints output result to the console.

        '''
        # ERRORS CHEKERS:
        print 'Retrieved rows:', b.shape[0], ' cols:', b.shape[1]
        if b.shape[1] != numcols or b.shape[0] != numrows:
            print 'WARNING:CALCULATED SHAPE DIFFERENT FROM EXPECTED!!'
            print '----This may cause problems in showing the data---'
            print "Calculated rows num: %s retrieved from file: %s" % (numrows, b.shape[0])
            print "Calculated columns num: %s retrieved from file: %s" % (numcols, b.shape[1])
            print "Consider adding the adjust manually the value with sr parameter"
            print "To correct use e.g. a.flip(sr=%i, sc=%i)" % (b.shape[0] - numrows, b.shape[1] - numcols)
        # VERBOSE MODE
        if v:
            print '------------------'
            print '------------------'
            print 'ORIGIN:', job_x_origin_mm, job_y_origin_mm
            print 'stage_step:', self.parameters.stage_step
            print 'RANGE X', self.parameters.rangeX
            print 'RANGE Y', self.parameters.rangeY
            print self.name
            print "num row:", numrows
            print "num cols:", numcols
            print "cols per rows", numrows * numcols
            print 'SHAPE:', b.shape
            
    def checkID(self):
        plt.ioff()
        if self.sample_infos != None:
                img = self.sample_infos.get_image()
                fig2t = plt.figure()
                ax1t = fig2t.add_subplot(111)
                ax1t.imshow(img)
                imgname = self.sample_infos.name
                label = self.sample_infos.label
                ax1t.set_title('IMGNAME: %s,%s SCAN: %s' %(imgname,label,self.name))
                from matplotlib.widgets import Button
                
                plt.subplots_adjust(bottom=0.2)
                def update(status):
                        self.logQC['MatchingID'] = status
                        plt.close(fig2t)
                        

                class Index:
                    ind = 0
                    def passed(self, event):
                        update('Passed')
                        
                
                    def failed(self, event):
                        update('Failed')
                
                callback = Index()
                axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
                axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
                bnext = Button(axnext, 'PASSED', color='g',hovercolor='y')
                bnext.on_clicked(callback.passed)
                bprev = Button(axprev, 'FAILED',color='r',hovercolor='y')
                bprev.on_clicked(callback.failed)
                plt.show()
                return True
           
                
        else:
            return False
                
    def QualityControl(self,user = 'GM'):
        folderpath = self.path[:-len(self.name)]
        from os import path
        
        if folderpath == '':
            folderpath = getcwd()
        
        files = listdir(folderpath)
        if 'QC-Failed.txt' in files:
            print "*********WARNING*********"
            print "File faild Quality control"
            raw_input("Press any key to contiune")
            status = 'Failed'
            self.logQC['result'] = status
        elif 'QC-Passed.txt' not in files:
            plt.switch_backend("Qt4Agg")
            try:
                self.checkID() 
            except IndexError:
                print "No image found!"
               
            self.plot(cm='dist',QC=True)
            self.SNR(QC=True,cmap='d')
            try:
                self.Total(QC=True,cmap='d')
            except IOError:
                print "No total found!"
            self.mfilter()
            self.subtractplane(QC=True)       
            self.plot(QC=True)
            status = 'Passed'
            if 'Failed' in self.logQC.values():
                status = 'Failed'
            self.logQC['user'] = user
            self.logQC['result'] = status
            import json
            json.dump(self.logQC,open("QC-%s.txt" %status,"w" ))
            print "Quality Control: %s" %(status)
            plt.switch_backend("TkAgg")
        elif 'QC-Passed.txt' in files:
            status = 'Passed'
            self.QCresult = status
            print "Quality Control: %s" %(status)
            self.logQC['user'] = user
            self.logQC['result'] = status
            
    def _get_QualityControl(self,user = 'GM'):
        folderpath = self.path[:-len(self.name)]
        if folderpath == '':
            folderpath = getcwd()
        files = listdir(folderpath)
        if 'QC-Failed.txt' in files:    
            self.logQC['result'] = 'Failed'
        elif 'QC-Passed.txt' not in files:
            self.logQC['result'] = 'Not Performed'
        elif 'QC-Passed.txt' in files:
            self.logQC['result'] = 'Passed'
        
        
    def checkmissing(self, replacewith=np.nan, repair=True, snr=False, sr=0):
        '''
        This is an important function that should be run every time you process
        data from microprofilometer.
        The .tag file is a series of number going form 1 to 255 then starting
        again form 0 to 255. For every measurment a tag is recorded if a
        measurment is missing the corrispetive tag is also missing (e.g. if you
        see in the tag file 14, 16 measurment 15 is missing).
        Missing measurments cause columns shifts, because missing elements are
        not replaced by any value a bad value is appended at the end of the
        measurment by default from the firmware.

        This function check if any missing value is present and if any are present
        and repair argument is set to True add a value (default is np.nan) for
        every missing measurment in the right position, and delate the last value.


        ------------------
         Keyword arguments:
             replacewith -- replace missing value with (default is np.nan)
             repair -- add missing value (default is True)
             sr -- Supplementary row have to be used in case
             snr -- save a corrected SNR matrix as CorrectSNR.npy matrix

        -----------------
        Returns:
            It wirte the self.missingvaluesarray parameter.

        '''
 
        if not path.exists(self.path + '.tag'):
            print 'File TAG do not exist'
            self.log['Missing Value correction'] = 'Not performed, no TAG file.'
            return
        from itertools import cycle
        if replacewith == 'lensmax':
            replacewith = self.lens.LENS_MAXd

        g = np.fromfile(self.path + '.tag', dtype='uint8')
        b = g.reshape(-1, self.parameters.numcols)
        lst = range(1, 256, 1)
        lst.append(0)  # add a 0 because it starts from 1,2,... , 255, 0
        lst = range(1, 256, 1)
        lst.append(0)

        processedmatrix = []
        for j in b:
            processedrow = []
            pool = cycle(lst)
            for i in j:
                kk = pool.next()
                if kk == i:
                    processedrow.append(False)
                else:
                    processedrow.append(True)
                    pool.next()

            processedmatrix.append(processedrow)

        missingvalue = np.array(processedmatrix)
        self.log['Number of missing value'] = np.count_nonzero(missingvalue)
        print 'Number of missng values: %s' % (self.log['Number of missing value'])

        if np.count_nonzero(missingvalue) == 0:
            self.log['Missing Value correction'] = 'No missing values to correct'

        if np.count_nonzero(missingvalue) > 0 and repair:
            print 'Re inizializing the data all processing will be lost....'
            self.inizialize()
            print 'Now correcting the missing value...'
            repariedmatrix = []
            for i, j in zip(self.array, missingvalue):
                # the zip returns zipped lines of the arrays
                leng = np.count_nonzero(j)  # count missing values in the line
                # get the indices of the missing values
                indexes = np.nonzero(j)
                # get the indices of last n/2 values
                todelete = indexes[0][-leng / 2:]
                corr = np.delete(i, todelete)  # delete the last elements
                # get the indices where to insert
                toinsert = indexes[0][:-leng / 2]
                # the missin values
                for d in toinsert:
                    corr = np.insert(
                        corr, d, replacewith)  # insert replacewith in
                    # the position of the missing value

                # when a missing value is in the last position this may cause
                # some problems because it is added another missing value in the
                # last position. So the row must be corrected.
                if len(corr) < len(i):
                    corr = np.append(
                        corr, np.full(
                            len(i) - len(corr), replacewith))
                if len(corr) > len(i):
                    sl = (len(i) - len(corr))
                    corr = corr[:sl]

                repariedmatrix.append(corr)
            print 'Reforming the final matrix...'
            b = np.array(repariedmatrix)
            rowtoflip = range(1, self.parameters.numrows, 2)
            b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
            self.array = b
            self.log['Missing Value correction'] = 'Performed'
            self.log['Number of missing value'] = np.count_nonzero(
                missingvalue)
            self.missingvaluesarray = missingvalue
    
    def interpolate_missing(self,method='cubic'):
        from scipy import interpolate
        import gc
        x = np.arange(0, self.array.shape[1])
        y = np.arange(0, self.array.shape[0])
        #mask invalid values
        #array = np.ma.masked_invalid(self.array)
        xx, yy = np.meshgrid(x, y)
        #get only the valid values
        gc.collect()
        self.array = interpolate.griddata((xx[~self.array.mask], yy[~self.array.mask]), 
                                  self.array[~self.array.mask].ravel(),
                                  (xx, yy),
                                   method=method)
    def savelog(self):
   
        self._mk_results_folder()
    

        with open(r'Results\postprocessing_log.txt', 'w') as f:
            f.write('file name: %s \n' % (self.path))
            f.write(
                'Added row to calculated from header: %s \n' %
                (self.log['supplementary rows']))
            f.write(
                'Added columns to calculated from header: %s \n' %
                (self.log['supplementary columns']))
            f.write(
                'Index of the column removed to fix the misalignment: %s \n' %
                (self.log['removed col']))
            f.write(
                'Method used to retrieve the array: %s \n' %
                (self.log['method']))

    def plot(self, vminx=None, vmaxx=None, unit='mm', data='', absolute=False,
             mode=False, plot=True, rectangle=[],save=False,cartesian=True,
            dilution=False, title = True, cm = None, QC = False ,
            orientation = 'auto',show_ROIs=False):
        '''
        This is the main function to plot the results. It handle most of the
        expeption and try to plot all the information about Lens and other
        parameters. So it will automatically see if a lens is present to retrive
        all the information about it. If not it will plot only the data
        regarding the measurment.

        The colorbar range can be chosen typing the lower and upper limit in
        the brakets (e.g. plot(17,21) will plot a color bar ranging form 17 to
        21), if None is left it will use the working range if a lens istance is
        present.

        Press 'c' twice to capture the coordinates of the current selection.

        ------------------
         Keyword arguments:
             vminx -- the lower colorbar limit (default is None)
             vmaxx -- the upper colorbar limit (default is None)
             unit --  the unit of the colorbar (default is 'mm')
             data -- A 2d matrix can be passed to be plotted (default is '')
             absolute -- if you need to use the coordinate system of the micor
             profilometer stages you can pass True to this argument, this argument
             is handy when you have done a scan of a full sample and you want to
             select a ROI to be analyzed with a higher resoltuion (default is False)
             auto -- If seted to True set the colorbar to min-max (default is False)



        '''
        def getax():
            x = plt.get(ax, 'xlim')
            y = plt.get(ax, 'ylim')
            return x, y

        def press(event):
            import sys
            print('press', event.key)
            sys.stdout.flush()
            if event.key == 'c':
                x, y = getax()
                print 'X coordinates: %0.4f, %0.4f Y cooridnates: %0.4f,%0.4f' % (x[0], x[1], y[0], y[1])

        # set colorbar orientation
        if self.parameters.colorbar_label == 'Roughness (microns)':
            labelres = 'patch_size'
        else:
            labelres = 'stage_step'
            
        if orientation == 'auto':
            if self.array.shape[0] < self.array.shape[1]:
                orientation = 'horizontal'
            else:
                orientation = 'vertical'
        
        if unit == 'um' or unit == 'micron' or unit == 'microns':
            self.parameters.colorbar_label = '($\mu m$) Repeatability=%s $\mu m$' % (
                self.lens.Repeatability)

        if (plot == True) and (self.lens != None) :
            self.parameters.xlabel = '(mm) \n  %s=%i $\mu m$ laser spot size=%s $\mu m$' % (
                labelres, self.parameters.stage_step * 1000, self.lens.X_laser_Spot_Size)
        else:
            self.parameters.xlabel = '(mm)'

        if self.log[
                'Best plane correction'] != 'Not performed' and vmaxx is None and vminx is None:
            mode = 2
        #This stretch colorbar values to min-max
        if mode == 'min-max':
            vminx = self.array.min()
            vmaxx = self.array.max()
        #You can also stre
        if mode == 1 or mode == 2 or mode == 3:
            stdar, meanar = np.std(self.array), np.mean(self.array)
            vminx = meanar - mode * stdar
            vmaxx = meanar + mode * stdar

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if cm == None:
            cmap = jet# Greys
        if cm == 'dist':
            minl = self.lens.LENS_MINd
            maxl =self.lens.LENS_MAXd
            midpoint = ( maxl - minl )/2
             
            cmap = LinearSegmentedColormap.from_list('mycmap', 
                                                     [(0,
                                                       'red'),
                                                    (0.5, 'green'),
                                                    (1,
                                                     'red')])
            cmap.set_over(color='k')
            cmap.set_under(color='w')
            vminx = minl
            vmaxx = maxl
        cm_class_str = "<class 'matplotlib.colors.LinearSegmentedColormap'>"
        if str(type(cm)) == cm_class_str :
            cmap = cm
        
        cmap.set_bad('w', 1.)  # set bad values
        if data == '':
            data = self.array
        #this make the origin to be in the bottom right corner 
        if cartesian:
            data = np.flipud(data)
        
        if dilution is not False:
            data = data[::dilution,::dilution]
            
        mesh = ax.pcolormesh(data, cmap=cmap, )
        ax.axis('scaled')
        ax.set_xlim(0, self.array.shape[1])
        ax.set_ylim(0, self.array.shape[0])
        
        if cartesian == False:
            ax.invert_yaxis()

        if absolute:
            unit = 'mm'
            start, end = ax.get_xlim()

            xticks = ticker.FuncFormatter(
                lambda x,
                pos: '{0:g}'.format(
                    x *
                    self.parameters.stage_step +
                    float(
                        self.parameters.job_x_origin_mm) +
                    self.parameters.offset))
            yticks = ticker.FuncFormatter(
                lambda x,
                pos: '{0:g}'.format(
                    x *
                    self.parameters.stage_step +
                    float(
                        self.parameters.job_y_origin_mm)))
            ax.xaxis.set_ticks(
                np.arange(
                    start,
                    end,
                    1 /
                    self.parameters.stage_step))
            ax.yaxis.set_ticks(
                np.arange(
                    start,
                    end,
                    1 /
                    self.parameters.stage_step))

        else:
            xticks = ticker.FuncFormatter(
                lambda x, pos: '{0:g}'.format(
                    x * self.parameters.stage_step))
            yticks = xticks


        ax.xaxis.set_major_formatter(xticks)
        ax.yaxis.set_major_formatter(yticks)
        ax.set_xlabel(
            '(mm) stage step=%i $\mu m$ ' %
            (self.parameters.stage_step * 1000))
        ax.set_ylabel('(mm)')
        if title:
            ax.set_title('%s' % (self.name))
        if self.lens != None:
            if self.parameters.colorbar_label != '($mm$) Repeatability=%s $\mu m$' % (
                self.lens.Repeatability):
                    unit = 'micron'

        if vminx is not None or vmaxx is not None:

            if unit == 'micron' or unit == 'um':
                mesh.set_clim(vmin=vminx, vmax=vmaxx)
                zticks = ticker.FuncFormatter(
                    lambda x, pos: '{0:g}'.format(x * 1000))
                cbar = fig.colorbar(
                    mesh,format=zticks,orientation=orientation,extend='both')
                cbar.set_label('($\mu m$)')
                if self.lens is not None:
                    ax.set_xlabel(self.parameters.xlabel)
                    cbar.set_label(self.parameters.colorbar_label)
            else:
                mesh.set_clim(vmin=vminx, vmax=vmaxx)
                cbar = fig.colorbar(mesh,orientation=orientation,extend='both')
                cbar.set_label('(mm)')
                if self.lens is not None:
                    ax.set_xlabel(self.parameters.xlabel)
                    cbar.set_label(self.parameters.colorbar_label)

        elif unit == 'micron' or unit == 'um':
            zticks = ticker.FuncFormatter(
                lambda x, pos: '{0:g}'.format(x * 1000))
            cbar = fig.colorbar(mesh, format=zticks, orientation=orientation)
            cbar.set_label('($\mu m$)$')
            if self.lens is not None:
                ax.set_xlabel(self.parameters.xlabel)
                cbar.set_label(self.parameters.colorbar_label)

        elif self.lens is not None and mode == False:
            cbar = fig.colorbar(mesh, orientation=orientation)
            cbar.set_label('(mm)')
            cbar.set_label(self.parameters.colorbar_label)
            ax.set_xlabel(self.parameters.xlabel)
            if data is not  '':
                mesh.set_clim(
                    vmin=self.lens.LENS_MINd,
                    vmax=self.lens.LENS_MAXd)

        fig.canvas.mpl_connect('key_press_event', press)
        
        if show_ROIs:
            #This works on the indices ROIs not on the ROIS saved
            for roi in self.ROIs_indices.keys():
                y1,y2,x1,x2 = self.ROIs_indices[roi]
                width = x2 - x1
                height = y2 -y1
                if cartesian:
                    y,x = self.array.shape
                    y2 = y -y2
                rectangle.append([x1,y2,width,height])
                print x1,y2,width,height
                #ax.annotate(roi,(x1 + width/2., y2),ha='center')
                ax.text(x1 + width/2., y2,roi,horizontalalignment='center',
                        verticalalignment='top',color = 'white',
                        bbox = dict(facecolor='black'), family='serif',
                        weight = 'bold')
        
        if rectangle != []:
            from matplotlib.patches import Rectangle
            try:
                ax.add_patch(
                    Rectangle(
                        (rectangle[0],
                         rectangle[1]),
                        rectangle[2],
                        rectangle[3],
                        edgecolor = 'k',
                        linewidth = 3,
                        facecolor='none'))
                print "ADDED"
            except:
                for i in rectangle:
                    ax.add_patch(
                    Rectangle(
                        (i[0],
                         i[1]),
                        i[2],
                        i[3],
                        edgecolor = 'k',
                        linewidth = 3,
                        facecolor='none'))
                print "ADDED"
        

                        
        if QC:
            from matplotlib.widgets import Button
            plt.subplots_adjust(bottom=0.2)
            def update(status):
                if cm == 'dist':    
                    self.logQC['workingdistance'] = status
                    plt.close(fig)
                    
                else:
                    self.logQC['plot-quality'] = status                  
                    plt.close(fig)
                    

            class Index:
                ind = 0
                def passed(self, event):
                    update('Passed')
                    
            
                def failed(self, event):
                    update('Failed')
           
            
            callback = Index()
            axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
            axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
            bnext = Button(axnext, 'PASSED', color='g',hovercolor='y')
            bnext.on_clicked(callback.passed)
            bprev = Button(axprev, 'FAILED',color='r',hovercolor='y')
            bprev.on_clicked(callback.failed)

        if save: self._save_plot(fig,'Plot %s.png' %(self.name))
        if plot:
            plt.show()
        else:
            return plt
            
    def profile_plot(self, col=None,row=None, start =None, stop =None,vminx=None,
                     vmaxx=None, unit='mm', absolute=False,mode=False, 
                     plot=True,save=False, title = True):
        '''
        plot profile
        '''
        if col != None and row == None:
            data = self.array[:,col][start:stop]
        elif row !=None and col == None:
            data = self.array[row][start:stop]
        else:
            print "Select a row or a column"
        
        if unit == 'um' or unit == 'micron' or unit == 'microns':
            self.parameters.colorbar_label = '($\mu m$) Repeatability=%s $\mu m$' % (
                self.lens.Repeatability)

        if (plot == True) and (self.lens != None) :
            self.parameters.xlabel = '(mm) \n  %s=%i $\mu m$ laser spot size=%s $\mu m$' % (
                'stage_step', self.parameters.stage_step * 1000, self.lens.X_laser_Spot_Size)
        else:
            self.parameters.xlabel = '(mm)'

        if self.log[
                'Best plane correction'] != 'Not performed' and vmaxx is None and vminx is None:
            mode = 2
        #This stretch colorbar values to min-max
        if mode == 'min-max':
            vminx = data.min()
            vmaxx = data.max()
        #You can also stre
        if mode == 1 or mode == 2 or mode == 3:
            stdar, meanar = np.std(data), np.mean(data)
            vminx = meanar - mode * stdar
            vmaxx = meanar + mode * stdar

        fig = plt.figure()
        ax = fig.add_subplot(111)
  
            
        ax.plot(data)

        if absolute:
            unit = 'mm'
            start, end = ax.get_xlim()

            xticks = ticker.FuncFormatter(
                lambda x,
                pos: '{0:g}'.format(
                    x *
                    self.parameters.stage_step +
                    float(
                        self.parameters.job_x_origin_mm) +
                    self.parameters.offset))
         
            ax.xaxis.set_ticks(
                np.arange(
                    start,
                    end,
                    1 /
                    self.parameters.stage_step))
         

        else:
            xticks = ticker.FuncFormatter(
                lambda x, pos: '{0:g}'.format(
                    x * self.parameters.stage_step))
           


        ax.xaxis.set_major_formatter(xticks)
        ax.set_xlabel(
            '(mm) stage step=%i $\mu m$ ' %
            (self.parameters.stage_step * 1000))
        ax.set_ylabel('(mm)')
        if title:
            ax.set_title('%s' % (self.name))
           
        if save: self._save_plot(fig,'Profile %s.png' %(self.name))
        if plot:
            plt.show()
        else:
            return plt
            
    def tredplotmayavi(self, array='array', show=True,outline=True,
                       axes=True,colorbar=True,unit='mm', 
                       cbarorientation='vertical'):
        if array is 'array':
            array = self.array
        
        if unit == 'micron' or unit == 'microns':
            array = array*1000

      

        from mayavi import mlab
        mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
        
        if hasattr(array, 'mask'):
            surH = mlab.surf(array, warp_scale='auto', mask=array.mask)
        else:
            surH = mlab.surf(array, warp_scale='auto')
        if colorbar: mlab.colorbar(surH, orientation=cbarorientation,
                                   title='Heights (%s)' %(cbunit))
        if outline: mlab.outline()
        if axes: mlab.axes()
    
        if show:
            mlab.show()
        else:
            return mlab

    def tredplotmayavimesh(self, array='array', show=True,outline=True,
                           axes=True,colorbar=True,unit='mm',
                           cbarorientation='vertical', scalars = 'array',save=False):
        if array == 'array':
            array = np.fliplr(self.array.copy())
            
        
                
        X1, Y1 = np.meshgrid(
            np.arange(0, round(array.shape[1] * self.parameters.stage_step, 3),
                      self.parameters.stage_step),
            np.arange(0, round(array.shape[0] * self.parameters.stage_step, 3),
                      self.parameters.stage_step))
                      
        if unit == 'micron' or unit == 'microns':
            array = array*1000
            X1 = X1*1000
            Y1 = Y1*1000
        if scalars == 'array':
            scalars = array
        elif scalars == '1order':
            scalars = self.subtractplane(order=1,array = array)
        elif scalars == '2order':
            scalars = self.subtractplane(order=2,array = array)           
        elif scalars == 'total':
            figs, scalars = self.Total(plot=False, matrix= True)
        from mayavi import mlab
        if hasattr(array, 'mask'):
            if array.mask.any() != False:
                print X1.shape, Y1.shape, array.shape
                mesh = mlab.mesh(X1[:array.shape[0],:array.shape[1]],
                                 Y1[:array.shape[0],:array.shape[1]], array,
                                 scalars=scalars, mask=array.mask)
            else:
                print X1.shape, Y1.shape, array.shape
                mesh = mlab.mesh(X1[:array.shape[0],:array.shape[1]],
                                 Y1[:array.shape[0],:array.shape[1]],
                             array, scalars=scalars)
                
        else:
            print X1.shape, Y1.shape, array.shape
            mesh = mlab.mesh(X1[:array.shape[0],:array.shape[1]],
                             Y1[:array.shape[0],:array.shape[1]],
                             array, scalars=scalars)
        mlab.colorbar(mesh, title='Heights (%s)' %(unit),orientation=cbarorientation)

        if outline: mlab.outline()
        if axes: mlab.axes()
        if show: mlab.show()
        if save: mlab.savefig(filename='%s.png' %(self.name))
        else:
            return mlab

    def tredplot(self, arr=None, stride='auto', plot=True, 
                 title='3D Plot',zprop=True,cm ='cw'):
        from matplotlib import cm
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax5 = fig.gca(projection='3d')
        ax5.set_aspect('equal')
        ax5.set_xlabel('(mm)')
        ax5.set_ylabel('(mm)')
        ax5.set_zlabel('(mm)')
        if cm == 'cw':
            colors = cm.coolwarm
        else:
            colors = cm.autumn
            
        colors.set_bad('k')

        if arr is None:
            arr = self.array.copy()

        if stride == 'auto':

            if arr.size > 250000:
                stride = arr.shape[0] / 500
                print 'Change stride to:', stride
            else:
                stride = 1

        X1, Y1 = np.meshgrid(
            np.arange(0, round(arr.shape[1] * self.parameters.stage_step, 3),
                      self.parameters.stage_step),
            np.arange(0, round(arr.shape[0] * self.parameters.stage_step, 3),
                      self.parameters.stage_step))

        if hasattr(arr, 'mask'):
            arr[arr.mask] = np.nan

        zmax = np.nanmax(arr)
        zmin = np.nanmin(arr)

        surf = ax5.plot_surface(X1, Y1, np.fliplr(arr), rstride=stride,
                                cstride=stride, cmap=colors, linewidth=0,
                                antialiased=True)
        surf.set_clim(vmin=zmin, vmax=zmax)
        if zmax < 1:
            zticks = ticker.FuncFormatter(
                lambda x, pos: '{0:g}'.format(x * 1000))
            cb = fig.colorbar(surf, format=zticks)
            cb.set_label('Heights ($\mu m$)')
            ax5.zaxis.set_major_formatter(zticks)
            ax5.set_zlabel('($\mu m$)')
            
            
        else:
            cb = fig.colorbar(surf)
            cb.set_label('Heights (mm)')
        
        #For keeping the Z axis with the same dimension 
        if zprop:
            max_range = np.array(
                [X1.max() - X1.min(), Y1.max() - Y1.min(), zmax - zmin]).max() / 2.0
    
            mid_x = (X1.max() + X1.min()) * 0.5
            mid_y = (Y1.max() + Y1.min()) * 0.5
            mid_z = (zmax + zmin) * 0.5
            ax5.set_xlim(mid_x - max_range, mid_x + max_range)
            ax5.set_ylim(mid_y - max_range, mid_y + max_range)
            ax5.set_zlim(mid_z - max_range, mid_z + max_range)

       
        ax5.set_xticks(ax5.get_xticks()[1:])
        ax5.set_yticks(ax5.get_yticks()[1:])
        ax5.set_title(title)
        ax5.view_init(elev=45., azim=125)
        if plot:
            plt.show()
        else:
            return plt
    
    def print_quality_index(self):
        if self.sample_infos != None:
            res = self.parameters.stage_step
            #float(self.sample_infos.width)/res * 
        
    def delcol(self, snr='.dist', v=False, col2remove=-1):
        '''
        It takes the name of the output file form the profilometer. It returns
        an array with the data rearranged. Deleting the last number of the row,
        this can solve the issue of misalingment of the rows.

        col2remove: default is -1 can be change 0 for the first one



        '''
        self.inizialize()
        b = np.delete(self.array, np.s_[col2remove], axis=1)
        rowtoflip = range(1, self.parameters.numrows, 2)
        b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
        self.diagnose(b, self.parameters.numcols - 1,
                      self.parameters.numrows, v)  # Run diagnose errors
        self.array = b
        self.log['method'] = 'delcol (for deleting a column)'
        self.log['removed col'] = col2remove

    def fliprows(self, snr='.dist', v=False):
        '''
        It takes the name of the output file form the profilometer. It returns
        an array with the data rearranged.


        '''
        rowtoflip = range(1, self.parameters.numrows, 2)
        b = self.array
        b[rowtoflip, :] = np.fliplr(
            self.array[rowtoflip, :])  # sobstitute rows
        self.diagnose(self.array, self.parameters.numcols,
                      self.parameters.numrows, v)  # Run diagnose errors
        self.array = b
        self.log['method'] = 'inverting the lines'
    
    def corners_detector(self, data = 'mask',w = 100,h = 100, 
                     Kernel_1 = (15,15),
                     Kernel_2 = (200,200),
                     pad_width = 100,
                     subdivide = True,
                     method = "cv2.TM_CCORR_NORMED",
                     plot = True
                     ):
        import cv2
        #Template
        if data == 'mask':
            imgo = self.array.mask
        else:
            imgo = data
            
        if h %2 != 0 or w%2 != 0:
            raise Warning('Template widht and height should be even.')
        template = np.ones((w,h),dtype=np.uint8)
        t_TR = template.copy()
        t_TL = template.copy()
        t_BR = template.copy()
        t_BL = template.copy()
        t_TL[h/2:,w/2:] = 0 #Top right
        t_TR[h/2:,:w/2] = 0 #Top left
        t_BL[:h/2,w/2:] = 0 # Bottom right
        t_BR[:h/2,:w/2] = 0 #Bottom left
        
        templates = {"TR":t_TR, "TL":t_TL,"BR": t_BR,"BL": t_BL}
        
        #Image preprocessing
        kernel = np.ones((Kernel_1[0],Kernel_1[1]),np.uint8)
        imgo = imgo.astype(np.float32)
        opening3 = cv2.morphologyEx(imgo, cv2.MORPH_OPEN, kernel)
        kernel2 = np.ones((Kernel_2[0],Kernel_2[1]),np.uint8)
        imgo2 = cv2.morphologyEx(opening3, cv2.MORPH_CLOSE, kernel2)
        imgop = np.pad(imgo2,pad_width=pad_width,mode='constant', constant_values=1)
        img2 = imgop.astype(np.uint8)
        corners =[]
        for j in templates:
            img = img2.copy()
            hz, wz = img.shape
            hp = 0 #this are the subdivision pad
            wp = 0
            #we subdivide the surface in 4 areas and we search for the template 
            #only in the are wher it should be
            if subdivide:
                if j == 'TL':
                    img = img[:hz/2,:wz/2]
                if j == 'TR':
                    wp = wz/2
                    img = img[:hz/2,wz/2:]
                if j == 'BL':
                    hp = hz/2
                    img = img[hz/2:,:wz/2]
                if j == 'BR':
                    wp = wz/2
                    hp = hz/2
                    img = img[hz/2:,wz/2:]
            template = templates[j]
            # Apply template Matching
            method2 = eval(method)
            res = cv2.matchTemplate(img,template,method2)
            min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
        
            # If the method is TM_SQDIFF or TM_SQDIFF_NORMED, take minimum
            if method2 in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
                top_left = min_loc
            else:
                top_left = max_loc
    
            corner = (top_left[0] + w/2 -pad_width + wp,
                      top_left[1] + h/2 -pad_width + hp)
            self.corners[j] = corner
            corners.append(corner)
        

        #PLOTTING
        if plot:
            import gc
            fig = plt.figure()
            axT = fig.add_subplot(221)
            axT.imshow(self.array) #plot results
            for i in corners:
                axT.scatter(i[0],i[1])
                
            fig2 = plt.figure()
            ax = fig2.add_subplot(221)     
            ax.imshow(imgo)
            ax.scatter(corners[0][0], corners[0][1])
            ax.set_xlim(corners[0][0]-180,corners[0][0]+180)
            ax.set_ylim(corners[0][1]-180,corners[0][1]+180)
            ax.set_title('Detected Point')
            ax.set_xticks([]) 
            ax.set_yticks([])
            
            ax2 = fig2.add_subplot(222)
            ax2.imshow(imgo)
            ax2.scatter(corners[1][0], corners[1][1])
            ax2.set_xlim(corners[1][0]-180,corners[1][0]+180)
            ax2.set_ylim(corners[1][1]-180,corners[1][1]+180)
            ax2.set_title('Detected Point')
            ax2.set_xticks([]) 
            ax2.set_yticks([])
        
            ax3 = fig2.add_subplot(223)
            ax3.imshow(imgo)
            ax3.scatter(corners[2][0], corners[2][1])
            ax3.set_xlim(corners[2][0]-180,corners[2][0]+180)
            ax3.set_ylim(corners[2][1]-180,corners[2][1]+180)
            ax3.set_title('Detected Point')
            ax3.set_xticks([]) 
            ax3.set_yticks([])
            
            ax4 = fig2.add_subplot(224)
            ax4.imshow(imgo)
            ax4.scatter(corners[3][0], corners[3][1])
            ax4.set_xlim(corners[3][0]-180,corners[3][0]+180)
            ax4.set_ylim(corners[3][1]-180,corners[3][1]+180)
            ax4.set_title('Detected Point')
            ax4.set_xticks([]) 
            ax4.set_yticks([])
            fig2.suptitle("Method:%s Position: %s" %(method,j))
            plt.show()
            gc.collect()

        
    
    def centroid_cor(self):
        '''
        This function realling the centroids and create a new matrix. The probe
        acquire signal while is running. The real center of the area that has
        been analyzed by the laser is hance shifted toward the scanning direction
        of a distance that we call integration distance (see fig. a.).
        This function correct this shift using spline, creating a matrix
        Xc (cfr. fig. a) interpolating the value from the array for every line
        and recopmputing the values for an allingn matrix (cfr. fig. b).

        a) x x x x x x     b) x x x x x
          x x x x x x         x x x x x
           x x x x x x   ->   x x x x x
          x x x x x x         x x x x x
           x x x x x x        x x x x x
        '''
        from scipy.interpolate import splrep, splev
        Xc, Y1 = np.meshgrid(
            np.arange(0, self.array.shape[1] * self.parameters.stage_step,
                      self.parameters.stage_step),
            np.arange(0, self.array.shape[0] * self.parameters.stage_step,
                      self.parameters.stage_step)
        )

        X1 = np.arange(0, self.array.shape[1] * self.parameters.stage_step,
                       self.parameters.stage_step)

        CCDFreq = float(self.parameters.CCDfreq)
        integrationdistance__mm = (1 / CCDFreq) * \
            float(self.parameters.X_axis_vel)
        # define the x positions of the centroids fig. a
        Xc[::2] = Xc[::2] + integrationdistance__mm / 2.0
        # define the x positions of the cnetroids  fig. a
        Xc[1::2] = Xc[1::2] - integrationdistance__mm / 2.0
        newyar = []
        for i, j in zip(Xc, self.array):
            tck = splrep(i, j, s=0)
            newyar.append(splev(X1, tck, der=0))
        self.array = np.array(newyar)

        tck2 = splrep(Xc[0], self.array[0], s=0)
        newy = splev(X1, tck2, der=0)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, aspect='equal')
        ax2.set_xlim(-self.parameters.stage_step,
                     self.parameters.stage_step * 5)
        ax2.set_ylim(-self.parameters.stage_step,
                     self.parameters.stage_step * 5)
        ax2.set_title('Interpolation first raw')
        ax2.set_xlabel('mm')
        ax2.set_ylabel('mm')
        ax2.plot(X1, newy, 'o', label='interpolated')
        ax2.plot(Xc[0], self.array[0], 'o', color='r', label='original')
        ax2.legend()

    def subtractplane(
            self,
            plot=False,
            inverted=True,
            eliminate_mask=False,
            dil=10,
            array=None,
            order=1,
            zprop=False,
            refplane=True,
            QC=False,
            returncalculated_plane= False
            ):
        '''
        See https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6
        for further documentation.
        '''
        import scipy.linalg
        from mpl_toolkits.mplot3d import Axes3D
        overwrite=False
        if array is None:
            array = self.array
            overwrite = True

        # regular grid covering the domain of the data Xr, Yr are undersampled
        # plane  to be plot faster
        X, Y = np.meshgrid(
            np.arange(
                0, array.shape[1], 1), np.arange(
                0, array.shape[0], 1))

        if hasattr(array, 'mask'):
            if array.mask.size !=1:
                # The masked values have to be extracted manually processing masked array will not
                # avoid computation over masked values
                x = X[~array.mask].flatten()
                y = Y[~array.mask].flatten()
                z = self.array[~array.mask].flatten()

            else:
                x = X.flatten()
                y = Y.flatten()
                z = array.flatten()

        else:
            x = X.flatten()
            y = Y.flatten()
            z = array.flatten()
        
        if order == 1:
            # best-fit linear plane
            A = np.c_[x, y, np.ones(len(x))]
            C, _, _, _ = scipy.linalg.lstsq(A, z)    # coefficients          
            # evaluate it on grid
            Z = C[0] * X + C[1] * Y + C[2]
            
        if order == 2:
            xy = np.dstack((x,y))[0]
            A = np.c_[np.ones(len(x)), xy,np.prod(xy, axis=1), xy**2]
            C,_,_,_ = scipy.linalg.lstsq(A, z)
            #Z = np.dot(np.c_[np.ones(x.shape), x, y, x*y, x**2, y**2], C).reshape(array.shape)
            Z = C[4]*X**2. + C[5]*Y**2. + C[3]*X*Y + C[1]*X + C[2]*Y + C[0]
        
        # store them in the object            
        self.parameters.plane_subtraction_coef = C
        # plot points and fitted surface
        xin= np.arctan((float(Z[0,1])-float(Z[0,0]))/self.parameters.stage_step)/np.pi*180
        yin= np.arctan((float(Z[1,0])-float(Z[0,0]))/self.parameters.stage_step)/np.pi*180
                 
        if plot or QC:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(X[::, ::dil], Y[::, ::dil], Z[
                            ::, ::dil], rstride=dil, cstride=dil, alpha=0.2)
            
            plt.xlabel('X')
            plt.ylabel('Y')
            ax.set_zlabel('Z')
            ax.axis('equal')
            ax.axis('tight')
            ax.set_title("Inclination: X %.2f$^\circ$ - Y  %.2f$^\circ$" %(xin,yin) )
            
            if refplane:
                Zref = np.zeros(self.array.shape)+np.mean(self.array.flatten())
                ax.plot_surface(X[::, ::dil], Y[::, ::dil], Zref[
                            ::, ::dil], rstride=dil, cstride=dil,
                            color='r',alpha=0.2)
             #For keeping the Z axis with the same dimension 
            if zprop:
                zmax = max(z)
                zmin = min(z)
                max_range = np.array(
                    [x.max() - x.min(), y.max() - y.min(), zmax - zmin]).max() / 2.0
        
                mid_x = (x.max() + x.min()) * 0.5
                mid_y = (y.max() + y.min()) * 0.5
                mid_z = (zmax + zmin) * 0.5
                ax.set_xlim(mid_x - max_range, mid_x + max_range)
                ax.set_ylim(mid_y - max_range, mid_y + max_range)
                ax.set_zlim(mid_z - max_range, mid_z + max_range)
            if QC:
                from matplotlib.widgets import Button
                plt.subplots_adjust(bottom=0.2)
                def update(status):   
                        self.logQC['Inclination'] = status
                        plt.close(fig)
    
                class Index:
                    ind = 0
                    def passed(self, event):
                        update('Passed')
                        
                
                    def failed(self, event):
                        update('Failed')
               
                
                callback = Index()
                axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
                axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
                bnext = Button(axnext, 'PASSED', color='g',hovercolor='y')
                bnext.on_clicked(callback.passed)
                bprev = Button(axprev, 'FAILED',color='r',hovercolor='y')
                bprev.on_clicked(callback.failed)
            plt.show()
        # Subtract plane
        
        
        print r"X axis inclination  %s degrees" %(xin)
        print r"Y axis inclination  %s degrees" %(yin)
        
        if inverted:
            if hasattr(self.array,"mask"):
                result = np.ma.array(Z - array.data)
                result.mask = array.mask
                
            else:
                result = Z - array
        else:
            if hasattr(self.array,"mask"):
                result = np.ma.array(array.data - Z)
                result.mask = array.mask
            else:
                result = array - Z
                

        if overwrite:
            self.array = result
            self.log['Best plane correction'] = 'Performed'
            
        if returncalculated_plane:
            return Z
        else:
            return result

    def manual_subtractplane(self, plot=False, inverted=True, dil=4,
                             useLens = False):
        '''
        See https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6
        for further documentation.
        '''
        import scipy.linalg
        from mpl_toolkits.mplot3d import Axes3D

        # regular grid covering the domain of the data Xr, Yr are undersampled
        # plane  to be plot faster
        X, Y = np.meshgrid(np.arange(0, self.array.shape[1], 1),
                           np.arange(0, self.array.shape[0], 1))

        x = []
        y = []
        if useLens:
            vminx = self.lens.LENS_MINd
            vmaxx = self.lens.LENS_MAXd
        else:
            stdar, meanar = np.std(self.array), np.mean(self.array)
            vminx = meanar - stdar * 2
            vmaxx = meanar + stdar * 2
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        mesh = plt.pcolormesh(self.array,vmin= vminx, vmax = vmaxx)
        ax.invert_yaxis()

        def onclick(event):
            global ix, iy
            ix, iy = event.xdata, event.ydata
            print 'x = %d, y = %d' % (
                ix, iy)
            x.append(int(round(ix)))
            y.append(int(round(iy)))
            return ix, iy

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        z = [self.array[i] for i in zip(y, x)]

        # best-fit linear plane
        A = np.c_[x, y, np.ones(len(x))]
        C, _, _, _ = scipy.linalg.lstsq(A, z)    # coefficients

        # evaluate it on grid
        Z = C[0] * X + C[1] * Y + C[2]

        # plot points and fitted surface
        if plot:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(X[::, ::dil], Y[::, ::dil], Z[
                            ::, ::dil], rstride=dil, cstride=dil, alpha=0.2)
            plt.xlabel('X')
            plt.ylabel('Y')
            ax.set_zlabel('Z')
            ax.axis('equal')
            ax.axis('tight')
            plt.show()
        # Subtract plane

        if inverted:
            result = Z - self.array
        else:
            result = self.array - Z
        self.array = result
        self.log['Best plane correction'] = 'Done manually'

    def sigmafilter(self, sigmas=1):
        meanar = np.mean(self.array)
        stdar = np.std(self.array)
        arrmask = np.ma.masked_outside(
            self.array,
            meanar - sigmas * stdar,
            meanar + sigmas * stdar)
        self.array = arrmask

    def mfilter(
            self,
            minz=None,
            maxz=None,
            badto=False,
            add=0,
            snr=True,
            lens=True,
            total='optimal',
            snrs=512):
        '''
        This is the mask filter that can be used to avoid bad values to be take into
        calulation it affects function as subtract bestplane, and plot every
        numpy function used will not take in account the masked value.

        It can maks also value that have a SNR less then argument snrs (default
        is 500 as suggested by Optimet).
        
        snr can also take snr matrix as input otherwise it will try to use the
        SNR matrix in the folde. 
        '''
        # The firs think to do is masking the np.nans
        self.array = np.ma.masked_invalid(self.array)
        log = ""
        if lens:
            if minz is None:
                minz = self.lens.LENS_MINd - add
    
            if maxz is None:
                maxz = self.lens.LENS_MAXd + add
    
            if minz == '-sigma':
                minz = (np.mean(self.array) - np.std(self.array))
    
            if maxz == '+sigma':
                maxz = (np.mean(self.array) + np.std(self.array))
    
            if minz == '-2sigma':
                minz = (np.mean(self.array) - 2 * np.std(self.array))
    
            if maxz == '+2sigma':
                maxz = (np.mean(self.array) + 2 * np.std(self.array))
    
            if minz == '-3sigma':
                minz = (np.mean(self.array) - 3 * np.std(self.array))
    
            if maxz == '+3sigma':
                maxz = (np.mean(self.array) + 3 * np.std(self.array))
    
            arrmask = np.ma.masked_outside(self.array, minz - add, maxz + add)
            log+= ' Between %s mm and %s mm' % (round(minz, 3), round(maxz, 3))
        else:
            arrmask = self.array
            
        if snr is not None:
            if snr is True:
                a = np.fromfile(self.path + '.snr', dtype=np.uint16)
                b = a.reshape(-1, self.parameters.numcols)
            else:
                b = snr
            rowtoflip = range(1, self.parameters.numrows, 2)
            b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
            arrmask2 = np.ma.masked_less(b, snrs)
            arrmask.mask = np.ma.mask_or(arrmask.mask, arrmask2.mask)
            log+= ' SNR > %i ' % (snrs)
        if total != False:
            try:
                tot =  np.fromfile(self.path + '.total', dtype=np.uint16)
                b = tot.reshape(-1, self.parameters.numcols)
                b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])
                if total == 'aceptable':
                    totmask = np.ma.masked_outside(b, 900, 18000)
                    log+="Total between 900 and 18000"
                    #note D.Blanco et al. / Precision Engineering 42(2015) 42-53
                    #Said total reccomanded 2000-16000
                if total == 'optimal':
                    totmask = np.ma.masked_outside(b, 1200, 16000)
                    log+="Total between 1200 and 16000"
                arrmask.mask = np.ma.mask_or(arrmask.mask, totmask.mask)
            except IOError:
                print "No total file found!"
            
        if badto == 'nan':
            arrmask[arrmask.mask] = np.nan
        elif badto is None or badto == 'None':
            arrmask[arrmask.mask] = None
        elif badto == 'Max':
            arrmask[arrmask.mask] = self.lens.LENS_MAXd
        elif badto == 'zero':
            arrmask[arrmask.mask] = 0
        self.array = arrmask
        self.log['Bad values filter'] = log 
        return log

# Data analyzing and instrument performance evaluation

    def ADF(
            self,
            vminx=None,
            vmaxx=None,
            auto=True,
            vlines=True,
            shift=0,
            plot=True,
            title=True):
        '''
        ADF, Amplitude Distribution Function
        '''
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        if title:
            ax2.set_title('%s' % (self.name))
        ax2.set_xlabel('Data (mm)')

        if vminx == '-sigma':
            vminx = (np.mean(self.array) - np.std(self.array))

        if vmaxx == '+sigma':
            vmaxx = (np.mean(self.array) + np.std(self.array))

        if vminx == '-2sigma':
            vminx = (np.mean(self.array) - 2 * np.std(self.array))

        if vmaxx == '+2sigma':
            vmaxx = (np.mean(self.array) + 2 * np.std(self.array))

        if vminx == '-3sigma':
            vminx = (np.mean(self.array) - 3 * np.std(self.array))

        if vmaxx == '+3sigma':
            vmaxx = (np.mean(self.array) + 3 * np.std(self.array))

        if self.log[
                'Best plane correction'] == 'Performed' and vmaxx is None and vminx is None:
            auto = False

        if vminx is not None or vmaxx is not None:
            auto = False
            bins = np.linspace(vminx, vmaxx, 100)
            if hasattr(self.array, 'mask'):
                ax2.hist(np.ravel(self.array[~self.array.mask]), bins=bins)
            else:
                ax2.hist(np.ravel(self.array), bins=bins)
            if vlines:
                ax2.vlines(
                    vminx - shift,
                    0,
                    ax2.axes.get_ylim()[1],
                    color='r',
                    linestyles=u'dashed')
                ax2.vlines(
                    vmaxx + shift,
                    0,
                    ax2.axes.get_ylim()[1],
                    color='r',
                    linestyles=u'dashed')
            if plot:
                plt.show()
            else:
                return plt

        if self.lens is not None and auto:
            bins = np.linspace(
                self.lens.LENS_MINd - 1,
                self.lens.LENS_MAXd + 1,
                1000)
            if hasattr(self.array, 'mask'):
                 x = ax2.hist(np.ravel(self.array[~self.array.mask]), bins=bins)
            else:
                 x = ax2.hist(np.ravel(self.array), bins=bins)
            print min(x[1]), max(x[1])
            print 'STDV: %s' % (np.std(np.ravel(self.array)))

            ax2.vlines(
                self.lens.LENS_MINd,
                0,
                ax2.axes.get_ylim()[1],
                color='r',
                linestyles=u'dashed')
            ax2.vlines(
                self.lens.LENS_MAXd,
                0,
                ax2.axes.get_ylim()[1],
                color='r',
                linestyles=u'dashed')

            if plot:
                plt.show()
            else:
                return plt

        else:
            ax2.hist(np.ravel(self.array))
            if plot:
                plt.show()
            else:
                return plt

    def timing_diagram(self, xlim='auto', plot=True):
        '''
        test the trigger and CCD frequency and plot a timing diagram
        '''
        from squarewave import sqwav
        time = (self.parameters.rangeX) / float(self.parameters.X_axis_vel)
        Freq = self.parameters.numcols / time
        Freq2 = 1 / (self.parameters.stage_step /
                     float(self.parameters.X_axis_vel))  # using velocity
        print 'Freq. from column:', Freq, 'Frq.:', Freq2, 'CCD:', self.parameters.CCDfreq
        integrationdistance__mm = (
            1 / self.parameters.CCDfreq) * float(self.parameters.X_axis_vel)
        real_x_resolution = integrationdistance__mm * 1000 + self.lens.X_laser_Spot_Size
        print 'Integration distance: %s mm Signal averaged area: %s micron' % (integrationdistance__mm, real_x_resolution)
        if Freq2 > 3000:
            print "Trigger frequency too high!"
        names = ['CCD', 'Trigger', ]
        # put always the lowest frequency at the end
        freq = [self.parameters.CCDfreq, Freq2, ]
        Triggerduration = 1000 / Freq2 - 0.0005
        durations = ['half period', Triggerduration]
        offsets = [0, 0.0005]
        a = zip(freq, durations, offsets, names)
        y_offset = 0
        for i in a:
            sqwav(i[0], i[1], i[2], i[3], y_offset=y_offset, xlim=xlim)
            y_offset += 1.5
        plt.yticks(np.arange(0.5, y_offset, 1.5), names)
        plt.ylim(-0.2, y_offset)
        self.log['Calculated spot diameter'] = real_x_resolution
        if plot:
            plt.show()
        else:
            return plt

    def laserspot_plot(
            self,
            centroid=True,
            angle=0,
            laserdirection='s',
            shifterror=False,
            plot=True):
        '''
        This function gives you an insight on the resolution it shows you the
        dimension of the laser spot size its initial and final postions, and
        the space that the laser has been moving while recording the line.
        The disalingment due to the CCD frequency can be evaluated.

        legend: True gives you a zoom on the first two measurment with the legend
        centroids: plot centroids in a separate plot
        '''
        from matplotlib.patches import Ellipse
        # Trigger meshgrid signlas
        X, Y = np.meshgrid(
            np.arange(0, self.parameters.stage_step * 5, self.parameters.stage_step),
            np.arange(0, self.parameters.stage_step * 5, self.parameters.stage_step)
        )
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_xlim(-self.parameters.stage_step,
                    self.parameters.stage_step * 5)
        ax.set_ylim(-self.parameters.stage_step,
                    self.parameters.stage_step * 5)
        ax.set_xlabel('mm')
        ax.set_ylabel('mm')
        height = self.lens.X_laser_Spot_Size / 1000.0
        width = self.lens.X_laser_Spot_Size / 500.0
        if laserdirection == 's':
            height = self.lens.X_laser_Spot_Size / 500.0
            width = self.lens.X_laser_Spot_Size / 1000.0
        CCDFreq = self.parameters.CCDfreq
        integrationdistance__mm = (1 / CCDFreq) * \
            float(self.parameters.X_axis_vel)
        # Final position signals
        X2, Y2 = np.meshgrid(
            np.arange(0, self.parameters.stage_step * 5, self.parameters.stage_step),
            np.arange(0, self.parameters.stage_step * 5, self.parameters.stage_step)
        )

        if shifterror:
            X2[::2] = X2[::2] + integrationdistance__mm
            X2[1::2] = X2[1::2] - integrationdistance__mm

        elif shifterror == False:
            X2[::2] = X2[::2] + integrationdistance__mm
            X[1::2] = X[1::2] + integrationdistance__mm

        XY = zip(X2.flatten(), Y2.flatten())
        # First plot of the laser spot
        for i in zip(X.flatten(), Y.flatten()):
            e = Ellipse(xy=i, width=width, height=height, alpha=0.7,
                        color='r', angle=angle)
            ax.add_patch(e)

        # This is the final position of the laser
        for x in XY:
            e = Ellipse(xy=x, width=width, height=height, alpha=0.4,
                        color='r', angle=angle)
            ax.add_patch(e)

        ax.plot(X, Y, "x", color='k')  # plot the trigger starting points
        # Explained plot
        figl = plt.figure()
        axl = figl.add_subplot(111, aspect='equal')
        axl.set_xlim(-self.parameters.stage_step,
                     self.parameters.stage_step * 5)
        axl.set_ylim(-self.parameters.stage_step,
                     self.parameters.stage_step * 5)
        axl.set_xlabel('mm')
        axl.set_ylabel('mm')
        axl.annotate(
            'Trigger start signal', xy=(0, self.parameters.stage_step * 4),
            xycoords='data',
            xytext=(-self.parameters.stage_step / 2,
                    self.parameters.stage_step * 4.35), textcoords='data',
            arrowprops={'arrowstyle': '->'})

        axl.annotate(
            'Laser spot', xy=(-width / 2, self.parameters.stage_step * 4),
            xycoords='data',
            xytext=(width / 2, self.parameters.stage_step * 4), textcoords='data',
            arrowprops={'arrowstyle': '<->'})

        axl.annotate(
            'Scanning step: %.3f mm' % (self.parameters.stage_step),
            xy=(0, self.parameters.stage_step * 4.25), xycoords='data',
            xytext=(self.parameters.stage_step, self.parameters.stage_step * 4.25),
            textcoords='data',
            arrowprops={'arrowstyle': '<->'})

        axl.annotate(
            'Exposure length: %.4f mm' % (integrationdistance__mm),
            xy=(0, self.parameters.stage_step * 3.75), xycoords='data',
            xytext=(integrationdistance__mm, self.parameters.stage_step * 3.75),
            textcoords='data',
            arrowprops={'arrowstyle': '<-'})

        axl.annotate(
            'Area in which the data are collected',
            xy=(-self.lens.X_laser_Spot_Size / 2.0, self.parameters.stage_step * 3.65),
            xycoords='data',
            xytext=(self.lens.X_laser_Spot_Size / 2000.0 + integrationdistance__mm,
                    self.parameters.stage_step * 3.65), textcoords='data',
            arrowprops={'arrowstyle': '<->'})

        axl.set_xlim(-self.lens.X_laser_Spot_Size / 2000.0 - 0.004, self.parameters.stage_step +
                     self.lens.X_laser_Spot_Size / 2000.0 + integrationdistance__mm + 0.0015)
        axl.set_ylim(0, self.parameters.stage_step * 6)
        for x in X[0]:
            for y in Y[:, 0]:
                e = Ellipse(xy=(x, y), width=width, height=height,
                            alpha=0.7, color='r', angle=angle)

                axl.add_patch(e)

        for x in XY:
            e = Ellipse(xy=x, width=width, height=height, alpha=0.4,
                        color='r', angle=angle)

            axl.add_patch(e)

        axl.plot(X, Y, "x", color='k')  # plot the trigger starting points

       # Centroid PLOT starts here:

        if centroid:

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111, aspect='equal')
            ax2.set_xlim(-self.parameters.stage_step,
                         self.parameters.stage_step * 5)
            ax2.set_ylim(-self.parameters.stage_step,
                         self.parameters.stage_step * 5)
            ax2.set_title('Centroids (o) and triggger signal (x)')
            ax2.set_xlabel('mm')
            ax2.set_ylabel('mm')
            X3, Y3 = np.meshgrid(np.arange(0, self.parameters.stage_step * 5,
                                           self.parameters.stage_step),
                                 np.arange(0, self.parameters.stage_step * 5,
                                           self.parameters.stage_step))
            if shifterror:
                X3[::2] = X3[::2] + integrationdistance__mm / 2.0
                X3[1::2] = X3[1::2] - integrationdistance__mm / 2.0

            elif shifterror == False:
                X3[::2] = X3[::2] + integrationdistance__mm / 2.0
                X3[1::2] = X3[1::2] + integrationdistance__mm / 2.0

            XY = zip(X3.flatten(), Y3.flatten())
            for x in XY:
                e = Ellipse(xy=x, width=integrationdistance__mm + width,
                            height=height, alpha=0.4, color='r', angle=angle)
                ax2.add_patch(e)

            # plot the trigger starting points
            ax2.plot(X, Y, "x", color='k', label='trigger signal')
            # plot the trigger starting points
            ax2.plot(X3, Y3, "o", color='r', label='centroids')
            plt.show()

    def ROIinteractive_select(self, scalingfactor=1,
                              useSNR=False,
                              save = 'ROI',
                              suffix='_ROI',
                              ):
        '''
        Return ROI of array to self.ROI.

        ------------------
         Keyword arguments:
             scalingfactor -- undesample matrix to allow visualization and avoid
             memory errors, notice the ROI will mantain the original size (default is 1)
             useSNR -- use SNR matrix to choose the ROI (default is False)



        -----------------
        Returns:
            missing values

        '''
        newarray = self.array
        counter = [1]
        from matplotlib.widgets import RectangleSelector
        if scalingfactor != 1:
            #from skimage.measure import block_reduce
            #newarray = block_reduce(self.array,(scalingfactor,scalingfactor))
            if useSNR:
                newarray = self.SNR(
                    matrix=True)[
                    ::scalingfactor,
                    ::scalingfactor]
                vmin, vmax = None, None
            else:
                newarray = self.array[::scalingfactor, ::scalingfactor]
                vmin, vmax = self.lens.LENS_MINd, self.lens.LENS_MAXd

        else:
            if useSNR:
                newarray = self.SNR(matrix=True)
                vmin, vmax = None, None
            else:
                newarray = self.array
                mean = np.mean(newarray)
                std = np.std(newarray)
                mode = 1
                vmin = mean - mode * std
                vmax = mean + mode * std
                #vmin, vmax = self.lens.LENS_MINd, self.lens.LENS_MAXd

        def line_select_callback(eclick, erelease):
            'eclick and erelease are the press and release events'
            x1, y1 = int(eclick.xdata), int(eclick.ydata)
            x2, y2 = int(erelease.xdata), int(erelease.ydata)
            print"(%s, %s) --> (%s, %s)" % (x1, y1, x2, y2)
            y1=y1*scalingfactor
            y2=y2*scalingfactor
            x1=x1*scalingfactor
            x2=x2* scalingfactor
            
            if save == 'ROI':
                self.ROI = self.crop(y1=y1,y2=y2,x1=x1,x2=x2,
                                     copy =True,
                                     suffix = suffix)
                print "Last ROI saved in self.ROI"
            if save == 'overwrite':
                self.crop(y1=y1,y2=y2,x1=x1,x2=x2,
                                     copy =False,
                                     )
            else: 
                self.ROIs_indices[suffix+str(counter[0])] = (y1,y2,x1,x2)
                counter[0] +=1
                print "All ROI indices saved in self.ROI_indices"
                                  
                                  

        def toggle_selector(event):
            print(' Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print(' RectangleSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print(' RectangleSelector activated.')
                toggle_selector.RS.set_active(True)

        print("\n      click  -->  release")

        fig, current_ax = plt.subplots()
        plt.imshow(newarray).set_clim([vmin, vmax])
        # drawtype is 'box' or 'line' or 'none'
        toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                               drawtype='box', useblit=True,
                                               # don't use middle button
                                               button=[1, 3],
                                               minspanx=5, minspany=5,
                                               spancoords='pixels',)
        plt.connect('key_press_event', toggle_selector)
        plt.show()

    
    def crop(self,y1,y2,x1,x2,copy=False,suffix=False):
                
        if copy:
            import copy
            surf = copy.deepcopy(self)
        else:
            surf = self
            
        surf.array = self.array[y1:y2,x1:x2]
        surf.parameters.numcols = y2-y1
        surf.parameters.numrows = x2-x1
        surf.parameters.rangeX = (x2-x1)*surf.parameters.stage_step
        surf.parameters.rangeY = (y2-y1)*surf.parameters.stage_step
        if suffix:
            surf.name = surf.name+suffix
            surf.special_ID = suffix
        
        if copy:
            return surf
        
    def crop_concentric_ROI(self,square_high__mm=1, usecorners = False,corner = 'TL'):
        #calculate the side of square in pts        
        pts = int(round(square_high__mm/self.parameters.stage_step/2,0))
        if usecorners:
                if self.corners[corner] != None: 
                    x1 = (self.corners['TL'][0] + self.corners['BL'][0])/2.
                    x2 = (self.corners['TR'][0] + self.corners['BR'][0])/2.
                    y1 = (self.corners['TL'][1] + self.corners['TR'][1])/2.
                    y2 = (self.corners['BL'][1] + self.corners['BR'][1])/2.
                    
                    xcent = int(round((x2 - x1)/2 + x1,0))
                    ycent = int(round((y2 - y1)/2 + y1,0))
 
                else: 
                    print "No corners found!"
                    return
            
        else:        
            ycent, xcent = np.array(self.array.shape)/2[0]
        
        self.array = self.array[ycent-pts:ycent+pts,
                                xcent-pts:xcent+pts]
        surf.parameters.numcols = y2-y1
        surf.parameters.numrows = x2-x1
        surf.parameters.rangeX = square_high__mm
        surf.parameters.rangeY = square_high__mm
        
   
            
    def concentric_ROI(self,step=20,plot='show',save=False):
        from scipy.stats.mstats import skew, kurtosis
        ycent, xcent = np.array(self.array.shape)/2
        distancefromcenter = []        
        Sq=[]
        Rskw=[]
        Rkurt=[]
        rectangles = []        
        
        i=1
        while step*i<xcent and step*i<ycent:
            ROI = self.array[ycent-step*i:ycent+step*i,xcent-step*i:xcent+step*i].flatten()
            print ROI.size
            distancefromcenter.append(step*i)
            Rskw.append(skew(ROI))
            Rkurt.append(kurtosis(ROI))     
            Sq.append(np.std(ROI)*1000)           
            rectangles.append((xcent-step*i,ycent-step*i,step*i*2,step*i*2))
            i+=1
            
        if plot !=False:
            distancefromcenter = np.array(distancefromcenter)*self.parameters.stage_step
            figx = plt.figure()
            axt = figx.add_subplot(111)
            axt.plot(distancefromcenter,Sq)
            axt.set_xlabel('Distance from the centered ROI (mm)')
            axt.set_ylabel('Sq ($\mu m$)')
            
            if plot == 'show':
                plt.show()         
            if save: self._save_plot(figx,'%s DistVsSq.png' %(self.name))
            
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.plot(distancefromcenter,Rskw,label= 'Skweness') 
            ax2.plot(distancefromcenter,Rkurt,label = 'Kurtosis')
            ax2.set_xlabel('Side of the centered ROI (mm)')
            ax2.set_ylabel('Normalized units')
            ax2.legend()
  
            if plot == 'show':
                plt.show()
                self.plot(rectangle=rectangles, unit='micron')
                
            if save:
                self._save_plot(fig2,'%s DistVsSkKur.png' %(self.name))
                ACD = self.plot(rectangle=rectangles, unit='micron',plot=False)
                ACD.savefig(r'%s\Results\%s ConcentricROIPlot.png'%(getcwd(),self.name))
                plt.close()
    
            
            
        return distancefromcenter,Sq,Rskw, Rkurt

    def margin_ROI(self,step=20,plot='show',save=False):
        from scipy.stats.mstats import skew, kurtosis
        import matplotlib.ticker as ticker 
        ycent, xcent = np.array(self.array.shape)/2
        margin = []        
        Sq=[]
        Rskw=[]
        Rkurt=[]
        rectangles = []        
        
        i=1
        while step*i<xcent and step*i<ycent:
            ROI = self.array[step*i:-step*i,step*i:-step*i].flatten()
            print ROI.size
            margin.append(step*i)
            Rskw.append(skew(ROI))
            Rkurt.append(kurtosis(ROI))     
            Sq.append(np.std(ROI)*1000)           
            rectangles.append((step*i,step*i,self.array.shape[1]-step*i*2,
                               self.array.shape[0]-step*i*2))
            print (step*i,step*i,self.array.shape[1]-step*i*2,
                               self.array.shape[0]-step*i*2)
            i+=1
            
        if plot !=False:
            margin = np.array(margin)*self.parameters.stage_step
            figx = plt.figure()
            axt = figx.add_subplot(111)
            axt.plot(margin,Sq)
            axt.set_xlabel('Margin (mm)')
            axt.set_ylabel('Sq ($\mu m$)')
            
            xticks = ticker.FuncFormatter(
                        lambda x,
                        pos: '{0:g}'.format(
                            x/self.parameters.stage_step))           
            
            axt2 = axt.twiny()
            axt2.set_xlabel('Margin (pt)')
            axt2.xaxis.set_major_formatter(xticks)
            
            if plot == 'show':
                plt.show()         
            if save: self._save_plot(figx,'%s DistVsSq.png' %(self.name))
            
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.plot(margin,Rskw,label= 'Skweness') 
            ax2.plot(margin,Rkurt,label = 'Kurtosis')
            ax2.set_xlabel('Margin (mm)')
            ax2.set_ylabel('Normalized units')
            axb2 = ax2.twiny()
            axb2.set_xlabel('Margin (pt)')
            axb2.xaxis.set_major_formatter(xticks)
            ax2.legend()
  
            if plot == 'show':
                plt.show()
                self.plot(rectangle=rectangles, unit='micron')
                
            if save:
                self._save_plot(fig2,'%s MarginVsSkKur.png' %(self.name))
                ACD = self.plot(rectangle=rectangles, unit='micron',plot=False)
                ACD.savefig(r'%s\Results\%s MarginROIPlot.png'%(getcwd(),self.name))
                plt.close()
    
            
            
        return margin,Sq,Rskw, Rkurt

    
    def analyize_ROI(self, bottomleftxy, width, height, show=False):
        '''

                     width
                    _______
                  h |     |
                    x,y___|
        '''
        from scipy.stats import skew
        ROI = self.array[
            xytopright[0]:xybottomleft[0],
            xytopright[1]:xybottomleft[1]]

        return np.std(ROI) * 1000, skew(ROI.flatten())

    def getROI(self, xytopleft, xybottomright,):
        '''
        Return ROI.
        '''
        ROI = self.array[
            xytopleft[1]:xybottomright[1],
            xytopleft[0]:xybottomright[0]]

        return ROI

    def calculate_profile_metrology(self,profile=None,roundx=3):
        ##  DEFINE SAMPLING LENGHT !!!!
    ####
        from scipy.stats.mstats import skew, kurtosis
        if self.log['Best plane correction']== 'Not Performed': 
            print 'WARINING BEST PLANE NOT PERFORMED!'
        if ROI == None:
            data=self.array
#        Amplitude parameters        
#    Rz, maximum height of the profile: defined on the sampling length: 
#        this parameter is frequently used to check whether the profile
#        has protruding peaks that might affect static or sliding contact
#        function.
#    
        Rz = round(np.max(profile)*1000,roundx)

#    Rv, maximum profile valley depth: depth of the deepest valley from 
#        the mean line, defined on the sampling length.
# 
        Rv = round(np.abs(np.min(profile))*1000,roundx)

#    Rt, total height of the profile: height between the deepest valley and 
#        the highest peak on the evaluation length.

        Rt = Rz + Rv

#    Ra, arithmetic mean deviation of the assessed profile: defined on 
#        the sampling length. Ra is used as a global evaluation of the roughness
#        amplitude on a profile. It does not say anything on the spatial
#        frequency of the irregularities or the shape of the profile.
#        Ra is meaningful for random surface roughness (stochastic) machined
#        with tools that do not leave marks on the surface, such as sand 
#        blasting, milling, polishing

        Ra = round(np.sum(np.abs(profile - np.mean(profile))/float(profile.size))*1000,roundx)
        
#    Rp, maximum profile peak height: height of the highest peak from the 
#        mean line, defined on the sampling length.
#    
        Rp = Rz - round(np.mean(profile)*1000,roundx)


#   
#    Rq, root mean square deviation of the assessed profile: corresponds to the
#        standard deviation of the height distribution, defined on the sampling
#        length. Rq provides the same information as Ra.
        Rq = round(np.std(profile)*1000,roundx)
#    
#    Rsk, skewness of the assessed profile: asymmetry of the height 
#        distribution, defined on the sampling length. This parameter is
#        important as it gives information on the morphology of the surface
#        texture. Positive values correspond to high peaks spread on a regular
#        surface (distribution skewed towards bottom) while negative values are
#        found on surfaces with pores and scratches (distribution skewed towards
#        top). It is therefore interesting when contact or lubrication functions
#        are required. However, this parameter does not give any information on 
#        the absolute height of the profile, contrary to Ra.
        
        Rsk = round(skew(profile.flatten()),roundx)
#    
#    Rku, kurtosis of the assessed profile: sharpness of the height 
#        distribution, defined on the sampling length.
#    
        Rku = round(kurtosis(profile.flatten()),roundx)
        
#    Rc, mean height of profile elements: defined on the evaluation length.
#        This parameter can be calculated on surfaces having texture cells or
#        grains. It is similar to the motif parameter R found in ISO 12085 and,
#        in that sense, it should be considered as a feature parameter 
#        (see ISO 25178).
        
        Rc = round(np.mean(profile)*1000,roundx)
        dicmetrology = {'Rz':Rz, 'Rv':Rv, 'Rt':Rt, 'Ra':Ra, 'Rp':Rp, 'Rq':Rq, 
                        'Rsk':Rsk, 'Rku':Rku, 'Rc':Rc}
        return dicmetrology
    

    
    def calculate_area_metrology(self,ROI=None,roundx=3,
                                 margin=None,subtractplane = False):
        from scipy.stats.mstats import skew, kurtosis
        if self.log['Best plane correction']== 'Not Performed': 
            print 'WARINING BEST PLANE NOT PERFORMED!'
        else:
            if subtractplane and self.log['Best plane correction']== 'Performed':
                print 'PLANE ALREADY SUBTRACTED!'
            else:
                print 'calculating...'
            
        if ROI == None:
            data=self.array
        if ROI != None:
            data=self.array[ROI[0]:ROI[1],ROI[2]:ROI[3]]
            
        if margin != None:
            data=data[margin:-margin,margin:-margin]
        if subtractplane:
            data = self.subtractplane(array = data)
#        Amplitude parameters        
#    Sz, maximum height of the profile: defined on the sampling length: 
#        this parameter is frequently used to check whether the profile
#        has protruding peaks that might affect static or sliding contact
#        function.
#    
        Sz = round(np.max(data)*1000,roundx)

#    Sv, maximum profile valley depth: depth of the deepest valley from 
#        the mean line, defined on the sampling length.
# 
        Sv = round(np.mean(data)-np.abs(np.min(data))*1000,roundx)



#    Sar, arithmetic mean deviation of the assessed profile: defined on 
#        the sampling length. Ra is used as a global evaluation of the roughness
#        amplitude on a profile. It does not say anything on the spatial
#        frequency of the irregularities or the shape of the profile.
#        Ra is meaningful for random surface roughness (stochastic) machined
#        with tools that do not leave marks on the surface, such as sand 
#        blasting, milling, polishing

        Sar = round(np.sum(np.abs(data - np.mean(data))/float(data.size))*1000,roundx)
        
#    Sp, maximum profile peak height: height of the highest peak from the 
#        mean line, defined on the sampling length.
#    
        Sp = Sz - round(np.mean(data)*1000,roundx)

#    St, total height of the profile: height between the deepest valley and 
#        the highest peak on the evaluation length.

        St = Sp + Sv


#   
#    Sq, root mean square deviation of the assessed profile: corresponds to the
#        standard deviation of the height distribution, defined on the sampling
#        length. Rq provides the same information as Ra.
        Sq = round(np.std(data)*1000,roundx)
#    
#    Ssk, skewness of the assessed profile: asymmetry of the height 
#        distribution, defined on the sampling length. This parameter is
#        important as it gives information on the morphology of the surface
#        texture. Positive values correspond to high peaks spread on a regular
#        surface (distribution skewed towards bottom) while negative values are
#        found on surfaces with pores and scratches (distribution skewed towards
#        top). It is therefore interesting when contact or lubrication functions
#        are required. However, this parameter does not give any information on 
#        the absolute height of the profile, contrary to Ra.
        
        Ssk = round(skew(data.flatten()),roundx)
#    
#    Sku, kurtosis of the assessed profile: sharpness of the height 
#        distribution, defined on the sampling length.
#    
        Sku = round(kurtosis(data.flatten()),roundx)
        
#    Sa, mean height of profile elements: defined on the evaluation length.
#        This parameter can be calculated on surfaces having texture cells or
#        grains. It is similar to the motif parameter R found in ISO 12085 and,
#        in that sense, it should be considered as a feature parameter 
#        (see ISO 25178).
        
        Sa = round(np.nanmean(np.abs(data))*1000,roundx)
        #Sa_peak = round(np.mean(data[data>0])*1000,roundx)
        
        dicmetrology = {'Sz_microns':Sz, 'Sv_microns':Sv, 'St_microns':St, 'Sar_microns':Sar, 'Sp_microns':Sp, 'Sq_microns':Sq, 
                        'Ssk':Ssk, 'Sku':Sku, 'Sa_microns':Sa}
        return dicmetrology
        
        

    def takenote(self):
        from scipy.stats import skew
        with open(r"C:\Users\OPdaTe\Documents\Microprofilometro - Acquisizioni\notes.txt", "a") as myfile:
            myfile.write(
                "%s -> %s microns , %s " %
                (self.name,
                 np.std(
                     self.array) *
                    1000),
                skew(
                    self.array.flatten()) +
                '\n')
    
    def go(self):
        self.mfilter()
        self.subtractplane()
        self.plot()
        
#♣Exporting tools.
    def GEOTIFF(self):
        import osgeo.gdal as gdal   
        """Array > Raster
        Save a raster from a C order array.

        :param array: ndarray

        +proj=tmerc +lat_0=51.4 +lon_0=7 +k=1 +x_0=0 +y_0=0 +ellps=WGS84
        +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
        """
        self._mk_results_folder()

        dst_filename = r'Results\ ' + r'%s.tiff' % (self.name)
        # You need to get those values like you did.
        x_min = float(self.parameters.job_x_origin_mm) + self.parameters.offset
        # x_min & y_max are like the "top left" corner.
        y_max = 300 - float(self.parameters.job_y_origin_mm)
        wkt_projection = 'WGS84'

        driver = gdal.GetDriverByName('GTiff')

        dataset = driver.Create(
            dst_filename,
            self.array.shape[1],
            self.array.shape[0],
            1,
            gdal.GDT_Float32, )

        dataset.SetGeoTransform((
            x_min,    # 0
            self.parameters.stage_step,  # size of the pixel...
            0,                      # 2
            y_max,    # 3
            0,                      # 4
            - self.parameters.stage_step))  # size of the pixel...

        dataset.SetProjection(wkt_projection)
        dataset.GetRasterBand(1).WriteArray(self.array)
        dataset.FlushCache()  # Write to disk.
        dataset = None  # Close dataset
        # return dataset, dataset.GetRasterBand(1)  #If you need to return,
        # remenber to return  also the dataset because the band don`t live
        # without dataset.

    def binary_export(self):
        '''
        Export of the flatten self.array parameter.
        '''

        self._mk_results_folder()

        with open(r'Results\%s.sdist' % (self.name), "wb") as f:
            newFileByteArray = bytearray(self.array.flatten())
            f.write(newFileByteArray)

    def npy_export(self):
        '''
        Export in .npy format, can be load with np.load()
        '''
        np.save(r'Results\%s' % (self.name), self.array)
        self.savelog()

    def mat_export(self, mdict = None, savemask =False):
        '''
        Export in .mat format.
        '''
        import scipy.io
 
        self._mk_results_folder()
        if mdict == None:
            if savemask:
                mdict = {'data' + self.name: self.array.data,
                         'mask' + self.name: self.array.mask}
            else:
                mdict = {'m' +
                    self.name: self.array}
        scipy.io.savemat(
            r'Results\ ' +
            self.name,
            mdict=mdict)
        self.savelog()

    def h5f_export(self, dataset=''):
        '''
        Export in .hdf5 format.
        '''
        import h5py
        if dataset == '':
            dataset = self.name
        h5f = h5py.File(r'Results\ '+dataset + '.hdf', 'w')
        h5f.create_dataset(dataset, data=self.array)
        self.savelog()
    
    def gwy_export(self,unit = "mm"):
        from gwyfile.objects import GwyContainer, GwyDataField, GwySIUnit
        obj = GwyContainer()
        obj["/0/data/title"] = self.name
        
        if  hasattr(self.array,"mask"):
            if unit == 'mm':
                data = GwyDataField(self.array.data.astype(np.float)/1000)
                data.si_unit_z = GwySIUnit(unitstr="m")
            if unit == 'micron' or unit == 'um':
                data = GwyDataField(self.array.data.astype(np.float)/1000)                
                data.si_unit_z = GwySIUnit(unitstr = "m")
            data.si_unit_xy = GwySIUnit(unitstr = "m")
            data.xres = self.array.shape[1]
            data.yres = self.array.shape[0]
            data.xreal = self.parameters.rangeX/1000
            data.yreal = self.parameters.rangeY/1000
            
            obj["/0/mask"] = GwyDataField(self.array.mask.astype(np.float))
            obj["/0/data"] = data
        else:
            if unit == 'mm':
                data = GwyDataField(self.array.astype(np.float)/1000)
                data.si_unit_z = GwySIUnit(unitstr = "mm")
            if unit == 'micron' or unit == 'um':
                data = GwyDataField(self.array.astype(np.float)/1000)                
                data.si_unit_z = GwySIUnit(unitstr = "m")
            data.xres = self.array.shape[1]
            data.yres = self.array.shape[0]
            data.xreal = self.parameters.rangeX/1000
            data.yreal = self.parameters.rangeY/1000
            data.si_unit_xy = GwySIUnit(unitstr = "m")
            obj["/0/data"] = GwyDataField(data)
        metadata = GwyContainer()
        metadata['OFFSET_mm'] = str(self.parameters.offset)
        metadata['RANGE_X_mm'] = str(self.parameters.rangeX)
        metadata['RANGE_Y_mm'] = str(self.parameters.rangeY)
        metadata['PROBE_FREQUENCY_hz'] = str(self.parameters.CCDfreq)
        metadata['JOB_X_ORIGIN_mm'] = str(self.parameters.job_x_origin_mm)
        metadata['JOB_Y_ORIGIN_mm'] = str(self.parameters.job_y_origin_mm)
        metadata['xAXIS_SPEED_mm/s'] = str(self.parameters.X_axis_vel)
        metadata['STAGE_STEP_mm'] = str(self.parameters.stage_step)
        metadata['RANGE_X_pt'] = str(self.parameters.numrows)
        metadata['RANGE_X_pt'] = str(self.parameters.numcols)
        metadata['LENS_NAME'] = self.lens.Name
        metadata['LENS_PN'] = self.lens.PN
        metadata['LASER_POWER'] = str(self.parameters.laserpower)
        metadata['LENS_SN'] = str(self.lens.SN)
        metadata['LENS_MINd'] = str(self.lens.LENS_MINd)
        metadata['LENS_MAXd'] = str(self.lens.LENS_MAXd)
        metadata['REPRODUCIBILITY'] = str(self.lens.Reproducibility)
        metadata['REPEATABILITY'] = str(self.lens.Repeatability)
        metadata['X_LASER_SPOT_SIZE'] = str(self.lens.X_laser_Spot_Size)
        obj["/0/meta"] = metadata
#        log = GwyContainer()
#        log[]
        obj.tofile("%s.gwy" %(self.name))
        
    def SPM_exportv2(self):

        self._mk_results_folder()
  
        header = ['ISO/TC 201 SPM data transfer format\n', 'general information\n',
                  '\n', '\n', '\n', '\n',
                  'Created by an image processing software.  Bogus acquisition parameters.\n',
                  'MAP_SC\n', '-1\n', '-1\n', '-1\n', '-1\n', '-1\n', '-1\n', '-1\n',
                  'scan information\n', 'REGULAR MAPPING\n',
                  'XYZ closed-loop scanner\n', 'sample XYZ scan\n', 'X\n',
                  'left to right\n', 'Y\n', 'top to bottom\n',
                  '151\n', '151\n', 'm\n', 'm\n', '0.0001\n', '0.0001\n',
                  'm\n', 'm\n', '0\n', '0\n', '0\n',
                  'm/s\n', '0.0\n', 'Hz\n', '0.0\n',
                  '\n', 'sample biased\n', '0.0\n', '0\n', '\n', '\n', '\n', '\n',
                  '\n', 'environment description\n', 'software\n', '300\n',
                  '1.0e5\n', '40\n', '\n', 'probe description\n', 'software\n',
                  '\n', '0.0\n', '0.0\n', '0.0\n', '0\n', '0\n', '0\n', '\n',
                  'sample description\n', 'Gray\n', '\n', '\n',
                  'single-channel mapping description\n', 'Gray\n', 'm\n', '\n',
                  'spectroscopy description\n', '\n', 'REGULAR\n',
                  '\n', 'n\n', '0.0\n', '0.0\n', '0.0\n', '0.0\n', '0\n', '0\n',
                  '\n', 'n\n', '0.0\n', '\n', 'data treatment description\n',
                  'post-treated data\n', '\n', '\n', '\n', '\n',
                  'multi-channel mapping description\n', '1\n', 'Gray\n',
                  'm\n', 'Gray\n', '\n', 'n\n', '\n', '\n', 'n\n', '\n', '\n',
                  'n\n', '\n', '\n', 'n\n', '\n', '\n', 'n\n', '\n', '\n', 'n\n',
                  '\n', '\n', 'n\n', '\n', '\n', 'n\n', '\n', '\n', 'n\n', '\n',
                  '\n', 'n\n', '\n', 'end of header\n']

        with open(r'Results\ ' + self.name + 'v2' + '.spm', 'w') as o:
            for i in header:
                o.write(i)

    def SPM_export(self):
        self._mk_results_folder()

        n = '\n'
        with open(r'Results\ ' + self.name + '.spm', 'w') as o:
            o.write('ISO/TC 201 SPM data transfer format' + n)
            o.write('general information' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write(
                'Created by Micropro UniVR.  Bogus acquisition parameters.' + n)
            o.write('MAP_SC' + n)
            o.write('-1' + n)
            o.write('-1' + n)
            o.write('-1' + n)
            o.write('-1' + n)
            o.write('-1' + n)
            o.write('-1' + n)
            o.write('-1' + n)
            o.write('scan information' + n)
            o.write('REGULAR MAPPING' + n)
            o.write('XYZ closed-loop scanner' + n)
            o.write('sample XYZ scan' + n)
            o.write('X' + n)
            o.write('left to right' + n)
            o.write('Y' + n)
            o.write('top to bottom' + n)
            o.write('%s' % (self.parameters.numrows) + n)
            o.write('%s' % (self.parameters.numcols) + n)
            o.write('m' + n)
            o.write('m' + n)
            o.write('%s' % (self.parameters.stage_step) + n)
            o.write('%s' % (self.parameters.stage_step) + n)
            o.write('m' + n)
            o.write('m' + n)
            o.write('0' + n)
            o.write('0' + n)
            o.write('0' + n)
            o.write('m/s' + n)
            o.write('%s' % (self.parameters.X_axis_vel) + n)
            o.write('Hz' + n)
            o.write('%s' % (self.parameters.CCDfreq) + n)
            o.write('' + n)
            o.write('sample biased' + n)
            o.write('0.0' + n)
            o.write('0' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('environment description' + n)
            o.write('software' + n)
            o.write('300' + n)
            o.write('1.0e5' + n)
            o.write('40' + n)
            o.write('' + n)
            o.write('probe description' + n)
            o.write('software' + n)
            o.write('' + n)
            o.write('0.0' + n)
            o.write('0.0' + n)
            o.write('0.0' + n)
            o.write('0' + n)
            o.write('0' + n)
            o.write('0' + n)
            o.write('' + n)
            o.write('sample description' + n)
            o.write('Gray' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('single-channel mapping description' + n)
            o.write('Gray' + n)
            o.write('m' + n)
            o.write('' + n)
            o.write('spectroscopy description' + n)
            o.write('' + n)
            o.write('REGULAR' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('0.0' + n)
            o.write('0.0' + n)
            o.write('0.0' + n)
            o.write('0.0' + n)
            o.write('0' + n)
            o.write('0' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('0.0' + n)
            o.write('' + n)
            o.write('data treatment description' + n)
            o.write('post-treated data' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('multi-channel mapping description' + n)
            o.write('1' + n)
            o.write('Gray' + n)
            o.write('m' + n)
            o.write('Gray' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('' + n)
            o.write('n' + n)
            o.write('' + n)
            o.write('end of header' + n)
            for j in self.array.flatten():
                o.write('%.7e\n' % (j))
    def pointCloud_generator(self,centroids=False,invertvalue=False):
        X1, Y1 = np.meshgrid(
        np.arange(0, round(self.array.shape[1] * self.parameters.stage_step, 3),
                  self.parameters.stage_step),
        np.arange(0, round(self.array.shape[0] * self.parameters.stage_step, 3),
                  self.parameters.stage_step))
        if centroids:
            CCDFreq = self.parameters.CCDfreq
            integrationdistance__mm = (
                1 / CCDFreq) * float(self.parameters.X_axis_vel)
            X1[::2] = X1[::2] + integrationdistance__mm / 2.0
            X1[1::2] = X1[1::2] - integrationdistance__mm / 2.0

        print 'X1 shape: ', X1.shape, 'Y1 shape: ', Y1.shape
        X = X1.flatten()
        print 'X1:', X.size
        Y = Y1.flatten()
        print 'y1:', Y.size
        Z = self.array.flatten()
        print 'z:', Z.size
        if invertvalue:
            Z = (self.lens.LENS_MAXd + 3) - Z
        xyz = zip(X, Y, Z)
        number_of_vertex = len(Z)
        return xyz, number_of_vertex       
        
    def ASCII_export(self, kind='ply', invertvalue=False, centroids=False):
        '''
        Method of ns class
        Export the array in a ASCII .ply format.

        Parameters
        ----------
        invertvalue : Boole
                 Peaks become valleys and viceversa. (default is False)
        lensfilter : Boole
                 Use the lens filter to avoid point out of range (NOT TESTED!) (default is False)
        centroids : Boole
                 Parameter used to alligned the centroid, if this taks is not perfomed automatically by the microprofilometer. (default is False)

        Returns
        -------
        Save a .ply file in the Result folder.
        '''
      
        n = '\n'
        xyz,nv = self.pointCloud_generator(centroids=centroids,invertvalue=invertvalue)
        self._mk_results_folder()
        if kind == 'ply':
            with open(r'Results\ ' + self.name + '.ply', 'w') as o:
                o.write('ply' + n)
                o.write('format ascii 1.0' + n)
                o.write('element vertex %s ' % (nv) + n)
                o.write('property float x' + n)
                o.write('property float y' + n)
                o.write('property float z' + n)
                o.write(
                    'comment Lens: %s SN: %s' %
                    (self.lens.Name, self.lens.SN) + n)
                o.write(
                    'comment stage_step: %s mm' %
                    (self.parameters.stage_step) + n)
                o.write(
                    'comment Repeatability: %s Laser Spot Size: %s ' %
                    (self.lens.Repeatability, self.lens.X_laser_Spot_Size) + n)
                o.write('end_header' + n)
                for i in xyz:
                    o.write('%s\t%s\t%s' % (i[0], i[1], i[2]) + n)

        if kind == 'xyz':
            with open(r'Results\ ' + self.name + '.XYZ', 'w') as o:
                o.write('# Channel: Gray' + n)
                o.write('# Lateral units: m' + n)
                o.write('# Value units: m' + n)
                for i in xyz:
                    o.write('%.4e\t%.4e\t%.4e' % (i[0], i[1], i[2]) + n)

    def exp_appyreport(self, metrology=3,flipyax=None):
        '''
        exp_appyreport export report with data directly from the project. The
        template used is called template.odt. The variables can be acessed
        using conditional fiedls Ctrl+F2 from the .odt document and placing 'true'
        in the 'Condition' field and the name of the variable in the 'Then' field.
        All the local variabels can be accessed e.g. self.name nested variables
        must be reassign e.g. stage_step = self.parameters.stage_step.

        '''
        import gc
        from appy.pod.renderer import Renderer
        from io import BytesIO

        print 'Plotting...'
        # Save reslts
        if self.log['Best plane correction'] == "Not performed":
            unit = "mm"
        else:
            unit = "micron"
        originalfigsize = plt.rcParams['figure.figsize']
        plt.rcParams['figure.figsize'] = (8,6)
        self.plot(mode=2, plot=False, unit=unit)
        figfile = BytesIO()
        try:
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
                
            plt.savefig(figfile, format='png', dpi=50)
            
        except MemoryError:
            print 'Size exceed memory limits,diluting data...'
            self.plot(mode=2, plot=False, unit=unit,dilution=2)
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
            figfile = BytesIO()
            plt.savefig(figfile, format='png', dpi=50)
            
        figfile.seek(0)  # rewind to beginning of file
        img = figfile.getvalue()  # extract string (stream of bytes))
        figfile.close()
        plt.close()
        gc.collect()
        print 'ADF...'
        # ADF Image
        self.ADF(
            '-3sigma',
            '+3sigma',
            shift=-
            np.std(
                self.array),
            plot=False)
        imghist = BytesIO()
        plt.savefig(imghist, format='png', dpi=50)
        imghist.seek(0)  # rewind to beginning of file
        hist = imghist.getvalue()  # extract string (stream of bytes))
        imghist.close()
        plt.close()
        print 'Attributes...'
        repeatability = self.lens.Repeatability
        stage_step = self.parameters.stage_step * 1000
        CCDfreq = self.parameters.CCDfreq/1000
        laserspotsize = self.lens.X_laser_Spot_Size
        lensType = self.lens.Name
        x_axis_speed = self.parameters.X_axis_vel
        laserpower = self.parameters.laserpower
        measurement_range = round(
            (self.lens.LENS_MAXd - self.lens.LENS_MINd), 2)

        if self.sample_infos is not None:
            description = self.sample_infos.description
            materials = self.sample_infos.materials
            Sample_name = self.sample_infos.name
            print 'Sample image..'
            plt.imshow(self.sample_infos.get_image())
            plt.tick_params(which='both', bottom='off', top='off',
                            labelbottom='off', labelleft='off',
                            left='off', right='off')
            figfiles = BytesIO()
            plt.savefig(figfiles, format='png', dpi=50)
            figfiles.seek(0)  # rewind to beginning of file
            img_smp = figfiles.getvalue()  # extract string (stream of bytes))
            figfiles.close()
            plt.close()

        else:
            Sample_name = self.name
        # SNR image
        print 'SNR...'
        self.SNR(plot=True)[0]
        figfile2 = BytesIO()
        try:
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
            plt.savefig(figfile2, format='png', dpi=50)
            
        except MemoryError:
            print 'Size exceed memory limits,diluting data...'
            self.SNR(plot=True,dilution=2)[0]
            if flipyax == -1:
                plt.gca().invert_yaxis()
            figfile2 = BytesIO()
            plt.savefig(figfile2, format='png', dpi=50)
        figfile2.seek(0)  # rewind to beginning of file
        img_SNR = figfile2.getvalue()  # extract string (stream of bytes))
        figfile2.close()
        plt.close()

        # Timing diagram
        self.timing_diagram(plot=False,)
        figfile3 = BytesIO()
        plt.savefig(figfile3, format='png', dpi=50)
        figfile3.seek(0)  # rewind to beginning of file
        img_TDG = figfile3.getvalue()  # extract string (stream of bytes))
        figfile3.close()
        plt.close()

        if metrology !=0:
            from scipy.stats.mstats import skew, kurtosis
            l = 125
            ycenter,xcenter = np.array(self.array.shape)/2
            while self.array.shape[0] < l*2*metrology or self.array.shape[1] < l*2*metrology:
                l-=20
            #dROI1 = self.array[400:550, 500:650]
            dROI1 = self.array[ycenter-l:ycenter+l, xcenter-l:xcenter+l]
            SqROI1 = round(np.std(dROI1) * 1000, 2)
            SkROI1 = round(skew(dROI1.flatten()), 2)
            KuROI1 = round(kurtosis(dROI1.flatten()), 2)
            figfile4 = BytesIO()
            gh = self.tredplot(dROI1, plot=0, title=1)
            gh.savefig(figfile4, format='png', dpi=60, bbox_inches='tight')
            figfile4.seek(0)  # rewind to beginning of file
            ROI1 = figfile4.getvalue()  # extract string (stream of bytes))
            figfile4.close()
            plt.close()

            if metrology >1:
                dROI2 = self.array[ycenter-l*2:ycenter-l, xcenter-l*2:xcenter-l]
                SqROI2 = round(np.std(dROI2) * 1000, 2)
                SkROI2 = round(skew(dROI2.flatten()), 2)
                KuROI2 = round(kurtosis(dROI2.flatten()), 2)
                figfile5 = BytesIO()
                gh = self.tredplot(dROI2, plot=0, title=2)
                gh.savefig(figfile5, format='png', dpi=60, bbox_inches='tight')
                figfile5.seek(0)  # rewind to beginning of file
                ROI2 = figfile5.getvalue()  # extract string (stream of bytes))
                figfile5.close()
                plt.close()

            if metrology >2:
                dROI3 = self.array[ycenter+l:ycenter+l*2, xcenter+l:xcenter+l*2]
                SqROI3 = round(np.std(dROI3) * 1000, 2)
                SkROI3 = round(skew(dROI3.flatten()), 2)
                KuROI3 = round(kurtosis(dROI3.flatten()), 2)
                figfile6 = BytesIO()
                gh = self.tredplot(dROI3, plot=0, title=3)
                gh.savefig(figfile6, format='png', dpi=60, bbox_inches='tight')
                figfile6.seek(0)  # rewind to beginning of file
                ROI3 = figfile6.getvalue()  # extract string (stream of bytes))
                figfile6.close()
                plt.close()
        
        else:
            metrology = 'nm'

        try:
            remove('%s.odt' % (self.path))
        except OSError:
            pass
         # Here you find the path of the template
        renderer = Renderer(r"C:\Users\OPdaTe\Documents\Microprofilometro - Acquisizioni\Template_new_metro%s.odt" %(metrology), locals(), '%s.odt' % (self.path))
        renderer.run()
        #restore original figsize
        plt.rcParams['figure.figsize'] =  originalfigsize 
        gc.collect()

    def viewer(self, absolute=False,complete=True):
        # Based on PythonQwt https://pypi.python.org/pypi/guiqwt

        sx = 0
        sy = 0
        if absolute:
            sx = self.parameters.job_x_origin_mm+self.parameters.offset
            sy = self.parameters.job_y_origin_mm
        plt.switch_backend("TkAgg")
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
        plot = win.get_plot()
        if complete:
            SNR = make.image(
                self.SNR(
                    plot=True,cartesian=False)[1],
                xdata=[
                    0 +
                    sx,
                    self.array.shape[1] *
                    self.parameters.stage_step +
                    sx],
                ydata=[
                    0 +
                    sy,
                    self.array.shape[0] *
                    self.parameters.stage_step +
                    sy],
                title="SNR",
                interpolation='nearest')
            
            plot.add_item(SNR, z=2)
            if self.missingvaluesarray is not None:
                rowtoflip = range(1, self.parameters.numrows, 2)
                b = self.missingvaluesarray
                b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
                missing = make.image(
                    b,
                    title="Missing",
                    xdata=[
                        0 +
                        sx,
                        self.array.shape[1] *
                        self.parameters.stage_step +
                        sx],
                    ydata=[
                        0 +
                        sy,
                        self.array.shape[0] *
                        self.parameters.stage_step +
                        sy],
                    interpolation="nearest")
                plot.add_item(missing, z=1)

        if hasattr(self.array, 'mask'):
            filteredarr = self.array.copy()
            filteredarr[self.array.mask] = np.nan
            filtered = make.image(
                filteredarr,
                title="Filtered",
                xdata=[
                    0 +
                    sx,
                    self.array.shape[1] *
                    self.parameters.stage_step +
                    sx],
                ydata=[
                    0 +
                    sy,
                    self.array.shape[0] *
                    self.parameters.stage_step +
                    sy],
                interpolation="nearest")
            plot.add_item(filtered, z=4)

        array = make.image(
            self.array,
            title="Data",
            xdata=[
                0 +
                sx,
                self.array.shape[1] *
                self.parameters.stage_step +
                sx],
            ydata=[
                0 +
                sy,
                self.array.shape[0] *
                self.parameters.stage_step +
                sy],
            interpolation="nearest")

        plot.add_item(array, z=3)
        win.exec_()
        plt.close('all')
        return win

    def bilateralFilter(self, d=6, sigmaColor=30, sigmaSpace=30):
        '''
        bilateralFilter(src, d, sigmaColor, sigmaSpace[, dst[, borderType]]) -> dst

        '''
        import cv2
        import copy
        if hasattr(self.array,"mask"):
            arr = copy.copy(self.array.data)
            print "masked data"
        else:
            arr = self.array
        filtered_surf = cv2.bilateralFilter(
            arr.astype(np.float32), d, sigmaColor, sigmaSpace)
        self.array = self.array - filtered_surf
        self.log['Processing'] = self.log['Processing'].join(
            'Bilateral filter d:%s, sigmaColor:%s, sigmaSpace:%s ' %
            (d, sigmaColor, sigmaSpace))
        return filtered_surf

    def MedianFilter(self, ksize):
        '''
        medianBlur(src, ksize[, dst]) -> dst
        '''
        import cv2
        filtered_surf = cv2.medianBlur(self.array, ksize)
        self.array = self.array - filtered_surf
        self.log['Processing'] = self.log['Processing'].join(
            'Median filter ksize:%s ' % (ksize))
        return filtered_surf

    def GaussianBlur(self, ksize, sigmaX):
        '''
        GaussianBlur(src, ksize, sigmaX[, dst[, sigmaY[, borderType]]]) -> dst
        '''
        import cv2
        filtered_surf = cv2.GaussianBlur(self.array, ksize, sigmaX)
        self.array = self.array - filtered_surf
        self.log['Processing'] = self.log['Processing'].join(
            'GaussianBlur ksize:%s, sigmaX: %s' % (ksize, sigmaX))
        return filtered_surf
    
    def normalize_0to255(self,astype = np.uint8):
        return (255*(self.array - np.min(self.array))/np.ptp(self.array)).astype(astype)

    def Segmentation(self, side, kind=None, mfilter=True,
                      overwrite=False,margin=None,norm=False):
        from scipy.stats import skew, kurtosis
        if kind is None:
            print " Select operation to perform on patches:"
            print " 1 - Roughenss"
            print " 2 - Skwenness"
            print " 3 - Kurtosis"
            print " 4 - False Color Composite"
            print " 5 - Equalized False Color Composite"
            kind = input('Type the number: ')

        
        if mfilter and hasattr(self.array, 'mask'):
            filteredarr = self.array.copy()
            filteredarr[self.array.mask] = np.nan
        else:
            filteredarr = self.array
            
        if margin != None:
            filteredarr=filteredarr[margin:-margin,margin:-margin]
        
        h, w = filteredarr.shape
        restor = h % side
        restoc = w % side
        
        if restoc != 0 and restor != 0:
            final = filteredarr[:-restor, :-restoc]
            print 'Eliminated %i columns and %i row' % (restoc, restor)
        elif restoc != 0:
            final = filteredarr[::, :-restoc]
            print 'Eliminated %i columns ' % (restoc)
        elif restor != 0:
            final = filteredarr[:-restor, ::]
            print 'Eliminated  %i row' % (restor)

        else:
            final = filteredarr

        numcolsa = final.shape[1] / side
        print numcolsa
        numrowsa = final.shape[0] / side
        print numrowsa
        
        #Form here starts the segmentation
        ini = 0
        arr = []
        arr2 = []
        R = []
        G = []
        B = []
        if kind == 1:
            if overwrite:
                self.parameters.colorbar_label = 'Roughness (microns)'
                self.name += "Seg"
            for i in range(1, numrowsa + 1):
                a = np.split(
                    np.ravel(
                        np.column_stack(
                            (final[ini:(side * i),::]))),
                    numcolsa)
                for k in a:
                    arr.append(np.nanstd(k))
                ini += side
           

        elif kind == 2:
            if overwrite:
                self.parameters.colorbar_label = 'Skwenness (microns)'
            for i in range(1, numrowsa + 1):
                a = np.split(
                    np.ravel(
                        np.column_stack(
                            (final[ini:(side * i),::]))),
                    numcolsa)
                for k in a:
                    arr.append(skew(k))
                ini += side

        elif kind == 3:
            if overwrite:
                self.parameters.colorbar_label = 'Kurtosis (microns)'
            for i in range(1, numrowsa + 1):
                a = np.split(
                    np.ravel(
                        np.column_stack(
                            (final[ini:(side * i),::]))),
                    numcolsa)
                for k in a:
                    arr.append(kurtosis(k))
                ini += side

        elif kind == 4 or kind == 5:
            for i in range(1, numrowsa + 1):
                a = np.split(
                    np.ravel(
                        np.column_stack(
                            (final[ ini:( side * i),::]))),
                    numcolsa)
                for k in a:
                    R.append(np.std(k))
                    G.append(skew(k))
                    B.append(kurtosis(k))
                ini += side

            R_arr = np.array(R).reshape(numrowsa, numcolsa) - np.min(R)
            G_arr = np.array(G).reshape(numrowsa, numcolsa) - np.min(G)
            B_arr = np.array(B).reshape(numrowsa, numcolsa) - np.min(B)

            if kind == 4:
                rgb_uint8 = (
                    np.dstack(
                        (R_arr /np.max(R_arr),
                            G_arr /np.max(G_arr),
                            B_arr /np.max(B_arr))) * 255.999).astype(np.uint8)
            if kind == 5:
                from skimage import exposure
                eq_Red = exposure.equalize_hist(R_arr / np.max(R_arr))
                eq_Blue = exposure.equalize_hist(B_arr / np.max(B_arr))
                rgb_uint8 = (
                    np.dstack(
                        (eq_Red,
                         G_arr /
                         np.max(G_arr),
                            eq_Blue)) *
                    255.999).astype(np.uint8)
            from PIL import Image
            img = Image.fromarray(rgb_uint8)
            self._mk_results_folder()
            img.save(
                r'Results\RGB composite %s pixel%s.tiff' %
                (self.name, side))
            plt.imshow(img, interpolation='nearest')
            plt.show()
            return
        
        elif kind == 6:
            if overwrite:
                self.parameters.colorbar_label = 'Bigerelle (microns)'
            for i in range(1, numrowsa + 1):
                a = np.split(
                    np.ravel(
                        np.column_stack(
                            (final[ini:(side * i),::]))),
                    numcolsa)
                for k in a:
                    arr.append(abs(np.nanmin(k)))
                    arr2.append(np.nanmax(k))
                ini += side
                
        elif kind == 7:
            #J. Europ. Opt. Soc. Rap. Public. 9, 14032 (2014)
            #Calculated in subareas
            #Divide area in subareas
        
            square_dist = np.zeros(numrowsa*numcolsa)
            subarea_mat=[]
            for i in range(1, numrowsa + 1):
                    a = np.split(
                        np.ravel(
                            np.column_stack(
                                (final[ini:(side * i),::]))),
                        numcolsa)
                    
                    for sub_ar in a:
                        subarea_mat.append(sub_ar.reshape(side,side))
                        
            x1 = random.randint(0,side-1)
            y1 = random.randint(0,side-1)
            x2 = random.randint(0,side-1)
            y2 = random.randint(0,side-1)
            
            dist_vector_r = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)*self.parameters.stage_step
            for n,ROI in enumerate(subarea_mat):            
                         square_dist[n] += math.sqrt(ROI[y1][x1] **2 + ROI[y1][x1] **2)

        
        if overwrite:
            self.parameters.stage_step = self.parameters.stage_step * side
            self.array = np.array(arr).reshape(numrowsa, numcolsa)
            self.log['Processing'] = self.log['Processing'].join(
                'Roughness segmentation side:%s' % (side))
            self.parameters.numcols = self.array.shape[1]
            self.parameters.numrows = self.array.shape[0]
            if norm:
                self.array = self.array/np.mean(arr)
                self.log['Processing'] = self.log['Processing'].join(
                ' divided by array mean')
                print 'Divinding by mean...'
                if norm =='to1':
                    minval = np.nanmin(self.array)
                    print 'Normalizing...'
                    self.array = (self.array-minval)/(np.nanmax(self.array)-minval)
                    self.log['Processing'] = self.log['Processing'].join(' norm to 1')
            
        else:
            if arr2 == []:
                return np.array(arr).reshape(numrowsa, numcolsa)
            else:
                return np.array(arr).reshape(numrowsa, numcolsa),np.array(arr2).reshape(numrowsa, numcolsa)

    def multi_segmentation(self,kind,step,
                           plot=True,
                           save = True,
                           xlim = None,
                           ylim = None,
                           norm = False
                           ):
        '''
        1 - Std. Deviation Roughness
        
        2 - Skewness 
        
        3 - Kurtosis
        
        6 - This is a 2D implementation of the 1D algorithm proposed in:
        The multi-scale roughness analyses and modeling of abrasion with the grit
        size effect on ground surfaces. 
        Wear 286-287 (2012) 124-135
        
        Maxence Bigerelle et al.
        '''
        multiseg = {}
        maxside = self.array.shape[0]/2
        if kind == 6:
            multiseg[round(1*self.parameters.stage_step*1000,1)] = [self.array,self.array]
        else:
            #add the array itself
            multiseg[round(1*self.parameters.stage_step*1000,1)] =  [self.array]     
        if  step == 'factors':
            steps = [i for i in range(2,maxside) if maxside%i == 0]
        else:
            steps = range(2,maxside,step)
        for t in steps:
            scale = round(t*self.parameters.stage_step*1000,1)
            multiseg[scale] = self.Segmentation(kind=kind,side=t, norm = norm)
          
        
        if plot:
            indices = np.array(sorted(multiseg.keys())) 
            if kind == 6:
                fig = plt.figure()
                ax = fig.add_subplot(221)      
                yminh = multiseg[indices[3]][0].flatten()*1000
                ymaxh = multiseg[indices[3]][1].flatten()*1000
                data = np.vstack([yminh, ymaxh]).T
                bins = np.linspace(0, np.nanmax(data), 30)
                ax.hist(data, bins, alpha=0.7, label=['-Ymin', 'Ymax'])
                ax.set_xlabel('Roughness amplitude (micron)')
                ax.set_title('Empirical distributions \n for evaluation \n lenght of %s micron'%(indices[3])) 
                ax.legend()
                
                ax = fig.add_subplot(222)
                indices = np.array(sorted(multiseg.keys()))
                yminh = multiseg[indices[5]][0].flatten()*1000
                ymaxh = multiseg[indices[5]][1].flatten()*1000
                data = np.vstack([yminh, ymaxh]).T
                bins = np.linspace(0, np.nanmax(data), 30)
                ax.hist(data, bins, alpha=0.7, label=['-Ymin', 'Ymax'])
                ax.set_xlabel('Roughness amplitude (micron)')
                ax.set_title('Empirical distributions \n for evaluation \n lenght of %s micron'%(indices[5])) 
                ax.legend()
                
                ax = fig.add_subplot(223)
                indices = np.array(sorted(multiseg.keys()))
                yminh = multiseg[indices[10]][0].flatten()*1000
                ymaxh = multiseg[indices[10]][1].flatten()*1000
                data = np.vstack([yminh, ymaxh]).T
                bins = np.linspace(0, np.nanmax(data), 30)
                ax.hist(data, bins, alpha=0.7, label=['-Ymin', 'Ymax'])
                ax.set_xlabel('Roughness amplitude (micron)')
                ax.set_title('Empirical distributions \n for evaluation \n lenght of %s micron'%(indices[10])) 
                ax.legend()
                
                ax = fig.add_subplot(224)
                indices = np.array(sorted(multiseg.keys()))
                yminh = multiseg[indices[15]][0].flatten()*1000
                ymaxh = multiseg[indices[15]][1].flatten()*1000
                data = np.vstack([yminh, ymaxh]).T
                bins = np.linspace(0, np.nanmax(data), 30)
                ax.hist(data, bins, alpha=0.7, label=['-Ymin', 'Ymax'])
                ax.set_xlabel('Roughness amplitude (micron)')
                ax.set_title('Empirical distributions \n for evaluation \n lenght of %s micron'%(indices[15])) 
                ax.legend()
                plt.tight_layout()
                
                
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(311)
                indices = np.array(sorted(multiseg.keys()))
                ymin = [np.nanmean(multiseg[i][0].flatten()*1000) for i in indices]
                ymax = [np.nanmean(multiseg[i][1].flatten()*1000) for i in indices]
                ax2.plot(indices,ymin,'x',color='b',label='-Ymin')
                ax2.set_xlabel('Evaluation length (micron)')
                ax2.set_ylabel('-Minimal\n roughness \n amplitude \n (micron)')
                ax3 = fig2.add_subplot(312)
                ax3.plot(indices,ymax,'x',color='r',label='Ymax')
                ax3.set_xlabel('Evaluation length (micron)')
                ax3.set_ylabel('maximal\n roughness\n  amplitude \n (micron)')
                ax4 = fig2.add_subplot(313)
                ax4.set_xlabel('Evaluation length (micron)')
                ax4.set_ylabel('maximal \n range \n roughness\n amplitude  (micron)')
                ax4.plot(indices,np.array(ymax)+np.array(ymin),'x',color='r',label='Ymax')
                plt.tight_layout()
                
  
                
        return multiseg
        
    def SLIC(
            self,
            minz=None,
            maxz=None,
            n_segments=250,
            compactness=0.1,
            enforce_connectivity=True):
        '''
        Method of ns class
        This is a wrapper for SLIC Segmentation from skimage. It allows a clever
        segmentation of the image where the superpixels don't have fixed areas
        but follows the countours.
        This skimage.segmentation.slic works only with
        value between -1 and 1 for this reason a preprocessing to normalize the
        data to one is needed. User can specifiy the range in which redefy the
        data.

        If calcRoughness argument is true the roughness of the superpixel is
        plotted.

        Parameters
        ----------
        minz : float
                 This parameter can be used to define a the minimum hight to be
                 used. (default is None)
        maxz : float
                 This parameter can be used to define the minumum hight to
                 be used. (default is None)
        n_segments : int
                 Number of segments. (default is 250)
        compactness : float
                 Compactness (see the documentation) it is set to 0.1 for gray
                 values images. (default is 0.1)
        enforce_connectivity : Boole
                 (see skimage.segemtation.slic for more information)
                 (default is True)

        Returns
        -------
        Return an array with the roughness calculated for every superpixel.

        '''
        from skimage.segmentation import slic
        from skimage.segmentation import mark_boundaries
        if minz is None or maxz is None:
                # the array must be with indixes between -1 and
                # 1 for this reason is normalized to min max
            lower_norm = self.array - np.min(self.array)
            norm2one = lower_norm / np.max(lower_norm)
        else:
            lower_norm = self.array - minz
            norm2one = lower_norm / (maxz - lower_norm)

        segments_slic = slic(norm2one, n_segments=n_segments,
                             compactness=compactness,
                             enforce_connectivity=enforce_connectivity)

        G = mark_boundaries(norm2one, segments_slic)
        plt.imshow(G)
        plt.show()

        calcRoughnessANS = raw_input('Do you want to calulate roughness? y/n ')
        if calcRoughnessANS.lower() == 'y':
            for i in np.unique(segments_slic):
                segment = np.ma.masked_where(segments_slic != i, self.array)
                Roughness = np.std(segment)
                self.array[segment.nonzero()] = Roughness

        self.log['Processing'] = self.log['Processing'].join(
            'SLIC minz=%s,maxz=%s,n_segments=%s, compactness=%s, enforce_connectivity=%s' %
            (minz, maxz, n_segments, compactness, enforce_connectivity))
        
        calccraksANS = raw_input('Do you want to calulate craks? y/n ')
       
        if calccraksANS.lower() == 'y':
            for i in np.unique(segments_slic):
                segment = np.ma.masked_where(segments_slic != i, self.array)
                mean = np.mean(segment)
                segment[segment > mean] = 0
                self.array[segment.nonzero()] = segment

        self.log['Processing'] = self.log['Processing'].join(
            'SLIC minz=%s,maxz=%s,n_segments=%s, compactness=%s, enforce_connectivity=%s' %
            (minz, maxz, n_segments, compactness, enforce_connectivity))
        return self.array

#Utilities 

    def _mk_results_folder(self):
         if not path.exists('Results'):
            makedirs('Results')
            
    def _save_plot(self,fig,name):
        self._mk_results_folder()
        dpi = 300
        while True:
            try:
                fig.savefig(r'Results\ '+name, dpi = dpi)
                break
            except:
                dpi /= 2
        plt.close(fig)


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
                    print finalpath 
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
                print "%s%s: %s  (%s ID:%s) " %(suff,i.name,label,name,idx)
                
            else:
                print " %s  has no ID!" %(i.name)
                status = False
                noIDs+=1
        if datasetID != []:
            
            self.ID = list(set(datasetID))
            print "Found %s dataset:" %(len(self.ID))
            try:
                import LabDBwrapper as Lab 
                self.dataset_infos = Lab.dataset(self.ID[0],r"C:\Users\OPdaTe\Documents\LabDatabase\SamplesDB")
            except ImportError:
                "No LabWrapper found"

            print self.ID
        if not status or LabelsNotInName:
            print "******WARNING!!!!********"
            print "Samples without label: %s" %(noIDs)
            print "Samples with label not in name: %s" %(LabelsNotInName)  
            if raw_input("Continue? y/n").upper() == "Y":
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
                print "NO samples infos!"
        return surfaces_index 
    
    def get_surface_by_label(self,label):
        surfaces_index = []
        for index,surf in enumerate(self.surfaces):
            if surf.sample_infos != None:
                if surf.sample_infos.label == label:
                    surfaces_index.append(index)
            else:
                print "NO samples infos!"
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
                 print "Sample: %s :" %(i[1])
                 indexs = self.get_surface_by_ID(i[0])
                 if indexs == []:
                     print ' *** MISSING ***'
                 else:
                     res = False
                     if ignoreQC:#Sample is counted as done even if QC fail
                         res=True
                     for j in indexs:
                         space = '          '+' '*len(i[1])
                         name = self.surfaces[j].name
                         QC = self.surfaces[j].logQC['result']
                         print "%s %s QC: %s" %(space,name,QC)
                         if QC == 'Passed':
                             res = True
                     if res:
                         done += 1
                        
            self.time_evaluation()
            status = done /float(self.dataset_infos.totnumber)
            print "STATUS:{0:.0f}%".format(status*100)
            timeleft = (self.totaltime/(int(status*100)))*int((100 -status*100))
            print "TIME LEFT: %s" %(timeleft)
                            
        else:
            print "No dataset infos!"
        
    def fliprows(self):
        for i in self.surfaces:
            print i.name
            i.fliprows()
            print "------------------"
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
        print 'Laser power difference: %s'%(laserpowerdiff)
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
            print 'Shapes do not match!'
        
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
        with open(r'Sq_comparison.csv', 'ab') as f:
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
            print "Select a row or a column"
        
     
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
            h,m,s = map(int,duration.split(':'))
            t0 += datetime.timedelta(hours = h,minutes = m,seconds = s)
        print t0
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
        fields = self.surfaces[0].calculate_area_metrology(margin=margin).keys()
        labels = []
        values = []
        class_maker = {}
        i=1
        if filename == None:
            filename = r'Metrology_complete.csv'
        with open(filename, 'ab') as f:
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
                        writer.writerow([surf.name,label[0],label[1]]+diz2.values())
                    elif split == 'filename':
                        if delimiter !=None:
                            dataclass = surf.name.split(delimiter)
                        else:
                            dataclass = surf.name
                            
                        writer.writerow([surf.name,dataclass[-2],label[-1]]+diz2.values())
                    else:
                        if uselabel:
                            writer.writerow([surf.name,label]+diz2.values())
                        else:
                            writer.writerow([surf.name]+diz2.values())
                    
                    if groupbylabel:
                        if surf.sample_infos.label not in class_maker:
                            class_maker[surf.sample_infos.label] = i
                            i+=1
                else:                   
                    if cutlabel == None and split == False:
                        labels.append(surf.name)
                        writer.writerow([surf.name]+diz2.values())
                    
                    elif split == 'filename':
                        if delimiter !=None:
                            dataclass = surf.name.split(delimiter)
                        else:
                            dataclass = surf.name
                        writer.writerow([surf.name,dataclass[-2],label[-1]]+diz2.values())
                        
                    elif cutlabel != None:
                        label = surf.name[cutlabel[0]:cutlabel[1]]
                        labels.append(label)
                        writer.writerow([surf.name,label]+diz2.values())
                    
                values.append(diz2.values())
   
        if plot:
            print "Plotting the data it will take some time plt.show() to show"
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
                
        rang = range(-ty, ty+2)
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
            print 'FLIPPING ROWS OF THE SAMPLES:'
            self.fliprows()
        
        log = True
        for j, i in enumerate(self.surfaces):
            print 'SURFACE %s of %s' % (j + 1, len(self.surfaces))
            print 'Name: %s' %(i.name)
            if checkmissing:
                print 'Checking for missing values...'
                i.checkmissing(replacewith=0)
                if log: self.log.append('Checked missed value')
                             
            if mfilter:
                print 'Filtering bad values...'
                i.mfilter(total=total)
                if log: self.log.append('Bad value filter')                
                
            if crop:
                print 'Cropping concentric the matrix...'
                i.crop(crop[0],crop[1],crop[2],crop[3])
                if log: self.log.append('Cropped. Side: y1 %s y2 %s x1 %s x2 %s' 
                %(crop[0],crop[1],crop[2],crop[3]))
                
            if cropconcentric !=None:
                print 'Cropping concentric the matrix...'
                i.crop_concentric_ROI(cropconcentric)
                if log: self.log.append('Square Cropped. Side: %s mm' %(cropconcentric*2))
                    
            if bilateral:
                print 'Using bilateral filter ...'
                i.bilateralFilter()
                if log: self.log.append('Bilateral filter')
                
            if subtractplane:
                print 'Subtracting plane...'
                i.subtractplane(order=order)
                if log: self.log.append('Plane subtraction')                
                
            if secmfilter != False:
                print 'Applying second filter...'
                i.mfilter('-%s' %(secmfilter),'+%s' %(secmfilter),snr=False)
                if log: self.log.append('Second bad value filter:-%s +%s' %(secmfilter))
                    
            if timingdiagram:
                print 'Timing diagram...'
                pltX = i.timing_diagram(plot=False)
                if saveres:
                    self._save_plot(pltX,'%s TD' %(i.name))
            if appyreport:
                print 'Exporting report ...'
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
        print 'Average Measurement Repeatability: %s microns' %(averagestd)
        print 'Roughness Repeatability: %s microns'%(RoughnessRepetibility)
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
            scipy.io.savemat(
                r'Results\ ' +
                self.name,
                mdict=mdict)
        self.save_log()
        
    def exp_csvReport(self):
        import csv
        with open('DatasetReport.csv', 'wb') as f:
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
            print ('Working on surface %s of %s:' %(n,num_surf))
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
                    data.si_unit_z = GwySIUnit(unitstr = "μm")
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
                    data.si_unit_z = GwySIUnit(unitstr = "μm")
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
                    print "Flipped: %s" %(i.name)
                
        # create the HDF5 file
        f = h5py.File('%s_%s.h5' % (self.name, kind), "w")

        # set the attributes that are specific for the file
        f.attrs['HDF5_Version'] = h5py.version.hdf5_version
        f.attrs['h5py_version'] = h5py.version.version
        diz = self.dict_surfaces_by_name() #sort the surface by name
        for i in sorted(diz, key= diz.get):
            print i.name
            corners_mat = np.array([None,None]) 
            if i.corners["TL"] != None and i.corners["BR"] != None:
                print "Inside corner mat"
                #create a croner matrix
                corners_mat = np.array([[i.corners["TL"],i.corners["TR"]],
                                       [i.corners["BL"],i.corners["BR"]]])
                if i.corners["Rotation"] != None:
                    corners_mat = np.rot90(corners_mat,rotation)
            
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
                if corners_mat != None:
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
            print '--------------'
            print "%s: %s of %s" %(i.name,count,len(self.surfaces))
            print '--------------'
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
            print '--------------'
            print k
            print '--------------'
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
            with open('Averageymin.csv', 'wb') as f:
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
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax2.legend(handles, labels,prop={'size':8},loc=4)
                
        handles, labels = ax3.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax3.legend(handles, labels, prop={'size':8},loc=4)
                
        handles, labels = ax4.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax4.legend(handles, labels, prop={'size':8},loc=4)
        
        handles, labels = ax2l.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax2l.legend(handles, labels,prop={'size':8},loc=4)
                
        handles, labels = ax3l.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax3l.legend(handles, labels, prop={'size':8},loc=4)
                
        handles, labels = ax4l.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
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
            with open('ymin.csv', 'wb') as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indicescsv[0])
                for i in ymincsv:
                    writer.writerow(i)
                
            with open('ymax.csv', 'wb') as f:
                import csv
                writer = csv.writer(f)
                writer.writerow(indicescsv[0])
                for i in ymaxcsv:
                    writer.writerow(i)
                    
            with open('maxrange.csv', 'wb') as f:
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
            print "%s: %s of %s" %(i.name,count,len(self.surfaces))
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
            with open('Averagestd.csv', 'wb') as f:
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
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        colors = [colormap(i) for i in np.linspace(0, 0.9,len(ax2.lines))]
        for i,j in enumerate(ax2.lines):
            j.set_color(colors[i])        
        
        ax2.legend(handles, labels,prop={'size':8},loc=4)
        #FOR THE LOG PLOT   
        #sort both labels and handles by labels        
        handles, labels = ax3.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
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
            with open('multistd.csv', 'wb') as f:
                    import csv
                    writer = csv.writer(f)
                    writer.writerow(indicescsv[0])
                    for i in stdcsv:
                        writer.writerow(i)
            with open('multistd_errors.csv', 'wb') as f:
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
                    print "appended" + roi.name +suffix
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
                print "appemded" + surf.name +suffix
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
        for ROIkey in self.surfaces[surface_index].ROIs_indices.keys():
            y1,y2,x1,x2 = self.surfaces[surface_index].ROIs_indices[ROIkey]
            
            if savetoROIsDS:
               self.crop(y1,y2,x1,x2,toROIdataset=True,suffix=ROIkey)
               
            if cropall:
                self.crop(y1,y2,x1,x2,suffix=ROIkey)
                
            
            print "done for %s" %(ROIkey)
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
        with open(name, 'wb') as f:
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
        with open('Segmentation_distribution_%s.csv' %(name), 'wb') as f, open('Dataset Seg %s.csv' %(name), 'wb') as h:
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
                print 'Wrong SHAPE found!!'
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
           for b in self.surfaces:
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
        print 'Plotting...'
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
            print 'Size exceed memory limits,diluting data...'
            plotF = self.plot(mode=2, plot=False, unit='micron',dilution=2,
                              title = title, orientation = orientation)
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
            self._save_plot(plotF,'%s SM' %(surf.name))
            
        plt.close()
        gc.collect()
        print 'ADF...'
        # ADF Image
        if ADF:
            ADFv = surf.ADF(plot=False,title = title, vminx = '-3sigma',
                            vmaxx = '+3sigma')
            self._save_plot(ADFv,'%s ADF' %(surf.name))
            plt.close()
        
 

        if surf.sample_infos is not None:

            print 'Sample image..'
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111)
            ax2.imshow(surf.sample_infos.get_image())
            ax2.tick_params(which='both', bottom='off', top='off',
                            labelbottom='off', labelleft='off',
                            left='off', right='off')
            self._save_plot(plt,'%s im' %(surf.name))
            plt.close()

        # SNR image
        print 'SNR...'
        if SNR:
            plotSNR = surf.SNR(plot=True, cmap = 'd',title = title,
                               orientation = orientation)[0]
            try:
                if flipyax == -1:
                    plt.gca().invert_yaxis()
                    plt.gca().invert_xaxis()
                self._save_plot(plotSNR,'%s SNR' %(surf.name))
                
            except MemoryError:
                print 'Size exceed memory limits,diluting data...'
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
                print 'Size exceed memory limits,diluting data...'
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
        with open('%s.csv' %(name), 'wb') as f:
            import csv
            writer = csv.writer(f)
            writer.writerow([firstcell]+header)
            for j in range(len(data)):
                writer.writerow([surfname[j]]+data[j])


#class super_dataset:    
#    def __init__(self):
#        self.datasets = []
#        
#        
#    def load_dataset(self):
#        folders = [i for i in listdir(getcwd()) if path.isdir(i) ]
#        
#    def load_mat_dataset(self,lateral_resolution,conversion_factor=1):
#        folders = [i for i in listdir(getcwd()) if path.isdir(i) ]
#        for j in folders:
#            surfd = surfacedataset()
#            import scipy.io
#            matfiles = [i for i in listdir(j) if i[-4:] == '.mat']
#            for i in matfiles:
#                surf = ns(None,init=False)
#                mat = scipy.io.loadmat('%s/%s' %(j,i))
#                name = 'm'+i.split('.')[0]
#                surf.array = np.array(mat[name])*conversion_factor
#                surf.name = name
#                surf.parameters.stage_step = lateral_resolution
#                surfd.surfaces.append(surf)
#            surfd.name = j
#            self.datasets.append(surfd)
#            
#    def multi_seg_stdv(self,step,save=False,average=True,savecsv=True):
#        indices = []
#        stddev = []
#        labels = []
#        for dataset in self.datasets:
#            ind, std = dataset.multi_seg_stddev(step,save=save,average=average)
#            indices.append(ind)
#            stddev.append(std)
#            labels.append(dataset.name)
#        figST = plt.figure()
#        axr = figST.add_subplot(111)
#        for x,y,l in zip(indices,stddev,labels):
#            axr.plot(x,y,label=l)
#            
#        if savecsv:
#            with open('multi_seg_stdv.csv', 'wb') as f:
#                import csv
#                writer = csv.writer(f)
#                for x,y,l in zip(indices,stddev,labels):                
#                    writer.writerow(list(l)+list(x))
#                    writer.writerow(list(l)+list(y))
#                    
#        
#        #To make sure we don't have twice the same color
#        colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
#        colorsx = [colormap(i) for i in np.linspace(0, 0.9, len(axr.lines))]
#
#        #linestyles = ['_', '-', '--', ':']
#        #import random
#        for t,j1 in enumerate(axr.lines):
#                j1.set_color(colorsx[t])
#        axr.set_xlabel('Evaluation lenght (micron)')
#        axr.set_ylabel('Roughness (microns)')
#        axr.legend()
#
#        return indices,stddev,labels

try:
    file2load = getcwd().split('\\')[-1]
    a = ns(file2load)
    a.fliprows()
    a.QualityControl()
#    a.viewer()
#    a.mfilter()
#    a.manual_subtractplane()
#    a.viewer(complete=0)
    print "----------------------INFOS---------------------------------------"
    print "Created surface called 'a' type: a. + TAB to see methods availabe."
    print "Type ?a.nameofthemethod() to get information about the method."
except:
    print "Failed to create surface from %s." % (file2load)
    print "Try a = ns('nameofthefile') for further details"
    try:
        d = surfacedataset()
        d.load_data_fromfolder()
        print "Created dataset called 'd' type d. + TAB to see methods availabe."
        del(a)
    except:
        print "Failed to create dataset!"
        del(d)

if __name__ == '__main__':
    a = ns('Gesso_new')
    a.fliprows()
    a.plot()
