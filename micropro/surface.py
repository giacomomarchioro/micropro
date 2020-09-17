"""
Spyder Editor
This is a script to elaborate data form ths profilometer.
"""
from __future__ import print_function
from __future__ import division 
if hasattr(__builtins__, 'raw_input'):
    input = raw_input
# Python 2.7 backward compatibility -end
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from os import listdir, getcwd, path, makedirs,remove, sep
from .lens import diclens
from .config import myconfig
__all__ = ['ns']

class ns:
    """New surface (ns) class.

    New surface class (ns):
    This class create an object of the micorprofilometer it allows different method
    to be used with the dot notation. The first thing to do is to create a new
    surface. This class automatically initialize the instances with a method.

    NOTE this module is design to scan over the X axis most of the feature works
    only with this set-up.
    Parameters
    ----------
    path : str
        the path of the .dist file.
    lens : lens obj, optional
        the lens used during the measurement, by default None
    init : bool, optional
        auto-initialization, by default True
    """

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
        self.name = path.split(sep)[-1]
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
                print('No files found!')


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
        """Loading acoustic microscopy data in mat file.

        Parameters
        ----------
        path : bool, optional
            Add the current working directory to the file., by default True

        Returns
        -------
        str
            Last string of the header
        """
        import scipy.io  as scio
        self.array = scio.loadmat(self.name)['data']
        base = self.name.split('_TOF.mat')[0].replace('-','_')
        self.array = np.ma.masked_invalid(self.array)
        if path:
            base = path.join(getcwd(),base)

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
        """
        This is a initializing function design to retrieve data form the header
        file of the microprofilometer. It returns all the parameters useful to
        perform further analysis of the data.
        A 2D array of raw data, un-flipped and unmasked is formed at the end
        of the process and stored in self.array. This will be the base for
        further processing.
        """
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

        print("WORKING ON:")
        print(self.name)
        #********************************************
        #**************HEADER VERSION 1**************
        #********************************************
        if 'RANGE_X' in self.header:
            import time
            print('Using VERSION 1 of the header')
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
                self.parameters.offset = eval(input(
                    'Please define offset manually: '))

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
                    lens = input(
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
                try:
                    from LabDBwrapper import Sample
                    p = myconfig['databasepath']
                    if p is not None and p != '':
                        self.sample_infos = Sample(int(self.header['JOB_ID']),p)
                except ImportError:
                    print("No labDB wrapper found.")
            self.parameters.finepower = float(
                self.header['PROBE_FINEPOWER'].replace(',', '.'))
            self.parameters.coursepower = float(
                self.header['PROBE_COURSEPOWER'].replace(',', '.'))


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
                try:
                    from LabDBwrapper import Sample
                    p = myconfig['databasepath']
                    if p is not None and p != '':
                        self.sample_infos = Sample(int(self.header['JOB_ID']),p)
                except ImportError:
                    print("No LabDB wrapper found")

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
                self.parameters.offset = eval(input(
                    'Please define offset manually: '))

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
                    lens = input(
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
                self.parameters.offset = eval(input(
                    'Please define offset manually: '))

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
        print(self.parameters.rangeX)
        self.parameters.numcols = int(
            self.parameters.rangeX / self.parameters.stage_step) + 1
        print('Num. rows: %s, Num. Cols: %s' % (self.parameters.numrows, self.parameters.numcols))

        # LOADING THE DATA
        arr = np.fromfile(self.path + '.dist', dtype=np.float32)
        print('Size expected: %s, Size retrieved: %s' % (self.parameters.numrows * self.parameters.numcols, arr.size))
        self.parameters.colorbar_label = '($mm$) Repeatability=%s $\mu m$' % (
            self.lens.Repeatability)

        if self.parameters.laserpower < 10:
            print("LASER POWER LOW, WAVELENGTH MIGHT BE UNSTABLE!")
        if self.parameters.numrows * self.parameters.numcols != arr.size:
            print("WARNING! SIZES DO NOT MATCH!")
            diff = self.parameters.numrows * self.parameters.numcols - arr.size
            print('Difference: %s' % (diff))
            print('Maybe rows: %s' % (diff / self.parameters.numrows))
            print('Maybe cols: %s' % (diff / self.parameters.numcols))
            try:
                b = arr.reshape(-1, self.parameters.numcols)
            except ValueError:
                g = eval(input('Add cols manually:'))
                b = arr.reshape(-1, self.parameters.numcols + g)
            srx = b.shape[0] - self.parameters.numrows
            scx = b.shape[1] - self.parameters.numcols
            print(" I am trying to fix it passing %s for rows and %s for columns" % (srx, scx))
            self.parameters.numrows += srx
            self.parameters.numcols += scx
            self.log['supplementary rows'] = srx
            self.log['supplementary columns'] = scx
            print("Now I am correcting the x-y range")
            self.parameters.rangeY += srx*self.parameters.stage_step
            self.parameters.rangeX += scx*self.parameters.stage_step

        else:
            b = arr.reshape(-1, self.parameters.numcols)
        self.array = b

        time = (self.parameters.rangeX) / float(self.parameters.X_axis_vel)
        Freq = self.parameters.numcols / time
        Freq2 = 1 / (self.parameters.stage_step /
                     float(self.parameters.X_axis_vel))  # using velocity
        print('Freq. from column:', Freq, 'Frq.:', Freq2, 'CCD:', self.parameters.CCDfreq)
        integrationdistance__mm = (
            1 / self.parameters.CCDfreq) * float(self.parameters.X_axis_vel)
        real_x_resolution = integrationdistance__mm * 1000 + self.lens.X_laser_Spot_Size
        self.log['Calculated spot diameter'] = real_x_resolution


    def Total(self, matrix=False, plot='show',
            cartesian=True, absolute=False, replacewith=0,dilution=False,
            QC=False,cmap='d',title=True, orientation = 'auto'
            ):
        """Compute the average Total value and allows displaying the Total matrix.

        Parameters
        ----------
        matrix : bool, optional
            reconstruct the total matrix otherwise it will compute the mean value, by default False
        plot : str, optional
            show the reconstructed Total matrix return a plot to show the data
            use the argument plot='show' or plt.show(), by default 'show'
        cartesian : bool, optional
            show the matrix using cartesian coordinate system 0,0 botttom left corner, by default True
        absolute : bool, optional
            use the absolute reference system of the stages, by default False
        replacewith : int, optional
            the value to be used instead of the missing values, by default 0
        dilution : bool, optional
            for saving memory the plot can be the result of a sampling, by default False
        QC : bool, optional
            perform the quality control, by default False
        cmap : str, optional
            the colormap to be used 'd' for showing discreate values 'c' for continuous, by default 'd'
        title : bool, optional
            the title to be displayed, by default True
        orientation : str, optional
            the orientation of the colorbar, by default 'auto'

        Returns
        ----------

        figs : matplotlib figures
            the created plot figures if plot is true
        b : the total matrix

        """
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
                print('Now correcting the missing value...')
                repariedmatrix = []
                for i, j in zip(b1, self.missingvaluesarray):
                    leng = np.count_nonzero(j)
                    indexes = np.nonzero(j)
                    todelete = indexes[0][int(round(-leng / 2)):]
                    corr = np.delete(i, todelete)
                    toinsert = indexes[0][:int(round(-leng / 2))]
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
                print('Reforming the final matrix...')
                b = np.array(repariedmatrix)
            else:
                b = a.reshape(-1, self.parameters.numcols)

            rowtoflip = list(range(1, self.parameters.numrows, 2))
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
                    cbar.set_ticklabels([900,'acceptable',1200,'optimal',16000,'acceptable',18000])




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
        print('Average Total: %s' % (self.log['Total']))
        return figs, b

    def SNR(self, matrix=False, plot='show',
            cartesian=True, absolute=False, replacewith=0,dilution=False,
            display='percentage',QC = False,cmap = 'c',title=True,
            orientation = 'auto' ):
        """Compute the average SNR value and allows displying the SNR matrix.

        Parameters
        ----------
        matrix : bool, optional
            if true reconstruct the SNR matrix, by default False
        plot : str, optional
            plot the reconstructed matrix, by default 'show'
        cartesian : bool, optional
            show the matrix using cartesian coordinate system 0,0 botttom left corner, by default True
        absolute : bool, optional
            use the absolute reference system of the stages, by default False
        replacewith : int, optional
            the value to be used instead of the missing values, by default 0
        dilution : bool, optional
            for saving memory the plot can be the result of a sampling, by default False
        display : str, optional
            shows the results as a percentage or as counts otherwise, by default 'percentage'
        QC : bool, optional
            perform the quality control, by default False
        cmap : str, optional
            the colormap to be used 'd' for showing discreate values 'c' for continuous, by default 'd'
        title : bool, optional
            the title to be displayed, by default True
        orientation : str, optional
            the orientation of the colorbar, by default 'auto'

        """

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
                print('Now correcting the missing values...')
                repariedmatrix = []
                for i, j in zip(b1, self.missingvaluesarray):
                    leng = np.count_nonzero(j)
                    indexes = np.nonzero(j)
                    todelete = indexes[0][-leng // 2:]
                    corr = np.delete(i, todelete)
                    toinsert = indexes[0][:-leng // 2]
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
                print('Reforming the final matrix...')
                b = np.array(repariedmatrix)
            else:
                b = a.reshape(-1, self.parameters.numcols)

            rowtoflip = list(range(1, self.parameters.numrows, 2))
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
                    tiks_lab = [30,'acceptable',50,'accurate',100]
                    labelKJ = 'SNR (%)'
                else:
                    mesh = axs.pcolormesh(b)
                    bounds = np.array([307,512,1024])
                    labelKJ = 'SNR (Raw value)'
                    tiks = [307,409,512,768,1024]
                    tiks_lab = [307,'acceptable',512,'accurate',1024]

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
        print('Average SNR: %s' % (self.log['SNR']))
        return figs, b

    def find_best_sphere(self):
        for line in self.array:
            pass


    def diagnose(self, b, numcols, numrows, v=False):
        """Diagnose the reconstructed array.

        Parameters
        ----------
        b : np.array
            The np.array reconstructed
        numcols : int
            The expected number of columns
        numrows : int
            The expected number of rows
        v : bool, optional
            If true a verbose diagnosis is used, by default False
        """
        # ERRORS CHEKERS:
        print('Retrieved rows:', b.shape[0], ' cols:', b.shape[1])
        if b.shape[1] != numcols or b.shape[0] != numrows:
            print('WARNING:CALCULATED SHAPE DIFFERENT FROM EXPECTED!!')
            print('----This may cause problems in showing the data---')
            print("Calculated rows num: %s retrieved from file: %s" % (numrows, b.shape[0]))
            print("Calculated columns num: %s retrieved from file: %s" % (numcols, b.shape[1]))
            print("Consider adding the adjust manually the value with sr parameter")
            print("To correct use e.g. a.flip(sr=%i, sc=%i)" % (b.shape[0] - numrows, b.shape[1] - numcols))
        # VERBOSE MODE
        if v:
            print('------------------')
            print('------------------')
            print('ORIGIN:', job_x_origin_mm, job_y_origin_mm)
            print('stage_step:', self.parameters.stage_step)
            print('RANGE X', self.parameters.rangeX)
            print('RANGE Y', self.parameters.rangeY)
            print(self.name)
            print("num row:", numrows)
            print("num cols:", numcols)
            print("cols per rows", numrows * numcols)
            print('SHAPE:', b.shape)

    def checkID(self):
        """Check if the ID used is the right one, confronting the photos in the database.

        This function is part of the quality control. It requires the user to compare the data in the
        database with the data collected.

        Returns
        -------
        bole
            True if there is an ID in the database (not necessarly the right one.)
        """
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
        """Perform a quality control on the measurement.

        This function allows to perform the quality control of the data collected. The user is asked to evaluate
        different parameters that are indicators of the quality of the scan.
        The results are stored in the self.logQC and in a separate file saved to disk QC-Passed.txt or QC-Failed.txt.

        Parameters
        ----------
        user : str, optional
            The name of the user who performed the quality control, by default 'GM'
        """
        folderpath = self.path[:-len(self.name)]
        from os import path

        if folderpath == '':
            folderpath = getcwd()

        files = listdir(folderpath)
        if 'QC-Failed.txt' in files:
            print("*********WARNING*********")
            print("File faild Quality control")
            input("Press any key to contiune")
            status = 'Failed'
            self.logQC['result'] = status
        elif 'QC-Passed.txt' not in files:
            plt.switch_backend("Qt4Agg")
            try:
                self.checkID()
            except IndexError:
                print("No image found!")

            self.plot(cm='dist',QC=True)
            self.SNR(QC=True,cmap='d')
            try:
                self.Total(QC=True,cmap='d')
            except IOError:
                print("No total found!")
            self.mfilter()
            self.subtractplane(QC=True)
            self.plot(QC=True)
            status = 'Passed'
            if 'Failed' in list(self.logQC.values()):
                status = 'Failed'
            self.logQC['user'] = user
            self.logQC['result'] = status
            import json
            json.dump(self.logQC,open("QC-%s.txt" %status,"w" ))
            print("Quality Control: %s" %(status))
            plt.switch_backend("TkAgg")
        elif 'QC-Passed.txt' in files:
            status = 'Passed'
            self.QCresult = status
            print("Quality Control: %s" %(status))
            self.logQC['user'] = user
            self.logQC['result'] = status

    def _get_QualityControl(self):
        """Get the results of the quality control.
        """
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


    def checkmissing(self, replacewith=np.nan, repair=True, snr=False):
        """Check if any valus is missing.

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

        Parameters
        ----------
        replacewith : values, optional
            The object to be used instead of the missing value, by default np.nan
        repair : bool, optional
            if Ture the matrix is corrected, by default True
        snr : bool, optional
            if Ture a corrected SNR matrix is saved as .npy array, by default False

        Returns
        -------
        It wirte the self.missingvaluesarray parameter.
        """
        if not path.exists(self.path + '.tag'):
            print('File TAG do not exist')
            self.log['Missing Value correction'] = 'Not performed, no TAG file.'
            return
        from itertools import cycle
        if replacewith == 'lensmax':
            replacewith = self.lens.LENS_MAXd

        g = np.fromfile(self.path + '.tag', dtype='uint8')
        b = g.reshape(-1, self.parameters.numcols)
        lst = list(range(1, 256, 1))
        lst.append(0)  # add a 0 because it starts from 1,2,... , 255, 0
        lst = list(range(1, 256, 1))
        lst.append(0)

        processedmatrix = []
        for j in b:
            processedrow = []
            pool = cycle(lst)
            for i in j:
                kk = next(pool)
                if kk == i:
                    processedrow.append(False)
                else:
                    processedrow.append(True)
                    next(pool)

            processedmatrix.append(processedrow)

        missingvalue = np.array(processedmatrix)
        self.log['Number of missing value'] = np.count_nonzero(missingvalue)
        print('Number of missng values: %s' % (self.log['Number of missing value']))

        if np.count_nonzero(missingvalue) == 0:
            self.log['Missing Value correction'] = 'No missing values to correct'

        if np.count_nonzero(missingvalue) > 0 and repair:
            print('Re inizializing the data all processing will be lost....')
            self.inizialize()
            print('Now correcting the missing value...')
            repariedmatrix = []
            for i, j in zip(self.array, missingvalue):
                # the zip returns zipped lines of the arrays
                leng = np.count_nonzero(j)  # count missing values in the line
                # get the indices of the missing values
                indexes = np.nonzero(j)
                # get the indices of last n/2 values
                todelete = indexes[0][int(round(-leng / 2)):]
                corr = np.delete(i, todelete)  # delete the last elements
                # get the indices where to insert
                toinsert = indexes[0][:int(round(-leng / 2))]
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
            print('Reforming the final matrix...')
            b = np.array(repariedmatrix)
            rowtoflip = list(range(1, self.parameters.numrows, 2))
            b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
            self.array = b
            self.log['Missing Value correction'] = 'Performed'
            self.log['Number of missing value'] = np.count_nonzero(
                missingvalue)
            self.missingvaluesarray = missingvalue

    def interpolate_missing(self,method='cubic'):
        """Interpolate the missing values with different methods.

        This method uses scipy interpolate, for computing the value where invalid values have been
        detected. Any masked values is replaced with the computed one self.mfilter() must be carried out
        before using this method. self.checkmissing() must be used before self.mfilter() for correcting
        any missing values.

        Parameters
        ----------
        method : str, optional
            The method used for interpolating the missing value, by default 'cubic'
        """
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
        """It saves the log in a .txt files.
        """

        self._mk_results_folder()


        with open(path.join("Results","postprocessing_log.txt"), 'w') as f:
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
        """plot the array of distances as color coded map.

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

        Parameters
        ----------
        vminx : float, optional
            the minimum value of the colorbar if None auto ranging is performed, by default None
        vmaxx : float, optional
            the maximum value of the colorbar if None auto ranging is performed by default None
        unit : str, optional
            the unit of the colorbar bar ('microns' or 'mm'), by default 'mm'
        data : np.array, optional
            the data to be plotted if '' the self.array is plotted, by default ''
        absolute : bool, optional
            if true the scanning reference system is used, by default False
        mode : bool, optional
            mode allows to strech the colorbar to one two or three sigmas (1,2,3) or to 'min-max', by default False
        plot : bool, optional
            when true the plot is distplayed, by default True
        rectangle : list, optional
            a list of vertices that can be used for overlaying rectangles (e.g. ROIs invistigated), by default []
        save : bool, optional
            if true save the resulting plot in the Results folder, by default False
        cartesian : bool, optional
            uses cartesian reference system (origin bottom left), by default True
        dilution : int, optional
            plot sampled values diluting the measurment by a factor, by default False
        title : str, optional
            the title of the plot, by default True
        cm : matplotlib.colormap, optional
            the color to be used in the plot if None the default colormap is used, by default None
        QC : bool, optional
            if True the setting for performing the quality control are used, by default False
        orientation : str, optional
            the oriantation of the colorbar ('vertical','horizontal') when auto the best oriantationis chosen
            depending on the aspect ratio, by default 'auto'
        show_ROIs : bool, optional
            shows the ROIs found as a set of rectangles., by default False

        Returns
        -------
        matplotlib.plt
            returns the matplotlib plot object
        """
        def getax():
            x = plt.get(ax, 'xlim')
            y = plt.get(ax, 'ylim')
            return x, y

        def press(event):
            import sys
            print(('press', event.key))
            sys.stdout.flush()
            if event.key == 'c':
                x, y = getax()
                print('X coordinates: %0.4f, %0.4f Y cooridnates: %0.4f,%0.4f' % (x[0], x[1], y[0], y[1]))

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
        cm_class_str = ["<class 'matplotlib.colors.LinearSegmentedColormap'>",
                        "<class 'matplotlib.colors.ListedColormap'>"]
        if str(type(cm)) in cm_class_str :
            cmap = cm
        else:
            cmap = plt.cm.viridis
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
            for roi in list(self.ROIs_indices.keys()):
                y1,y2,x1,x2 = self.ROIs_indices[roi]
                width = x2 - x1
                height = y2 -y1
                if cartesian:
                    y,x = self.array.shape
                    y2 = y -y2
                rectangle.append([x1,y2,width,height])
                print(x1,y2,width,height)
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
                print("ADDED")
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
                print("ADDED")



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
            return fig,ax

    def profile_plot(self, col=None,row=None, start =None, stop =None, unit='mm', absolute=False,mode=False,
                     plot=True,save=False, title = True):
        """plot a single profile form the distances array.

        Parameters
        ----------
        col : int, optional
            The coloumn to be ploted (if None a row must be selected), by default None
        row : int, optional
            The row to be plotted (if None a column must be selected), by default None
        start : int, optional
            The intial index of the column/row if None the start is the first index, by default None
        stop : int, optional
            The final index of the column/row if None the end is the last index, by default None
        unit : str, optional
            the unit of the colorbar bar ('microns' or 'mm'), by default 'mm'
        absolute : bool, optional
            if true the scanning reference system is used, by default False
        plot : bool, optional
            when true the plot is distplayed, by default True
        save : bool, optional
            if true save the resulting plot in the Results folder, by default False
        title : str, optional
            the title of the plot, by default True

        Returns
        -------
        matplotlib.plt
            the matplotlib plot object
        """
        if col != None and row == None:
            data = self.array[:,col][start:stop]
        elif row !=None and col == None:
            data = self.array[row][start:stop]
        else:
            print("Select a row or a column")

        if unit == 'um' or unit == 'micron' or unit == 'microns':
            self.parameters.colorbar_label = '($\mu m$) Repeatability=%s $\mu m$' % (
                self.lens.Repeatability)

        if (plot == True) and (self.lens != None) :
            self.parameters.xlabel = '(mm) \n  %s=%i $\mu m$ laser spot size=%s $\mu m$' % (
                'stage_step', self.parameters.stage_step * 1000, self.lens.X_laser_Spot_Size)
        else:
            self.parameters.xlabel = '(mm)'

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
        """3D plot of the distances array with enached distances z (not in scale).

        This function uses myavi library for producing a 3D plot of the distances array.

        Parameters
        ----------
        array : np.array, optional
            The array to be plotted using 'array' it plots self.array, by default 'array'
        show : bool, optional
            if True the plot is shown, by default True
        outline : bool, optional
            if True the outline of the plot is displyed, by default True
        axes : bool, optional
            if True the axes of the plot are shown, by default True
        colorbar : bool, optional
            if Ture the colorbar is shown, by default True
        unit : str, optional
            the unit of the colobar ('mm' or 'microns'), by default 'mm'
        cbarorientation : str, optional
            The orientation of the colorbar ('vertical' or 'horizontal'), by default 'vertical'

        Returns
        -------
        mlab
            mayavi mlab object
        """
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
        """3D plot of the distances array mantaining the proportions of the axes.

        This function uses myavi library for producing a 3D plot of the distances array.
        In this case the proportion of the three axes is mantained.

        Parameters
        ----------
        array : np.array, optional
            The array to be plotted using 'array' it plots self.array, by default 'array'
        show : bool, optional
            if True the plot is shown, by default True
        outline : bool, optional
            if True the outline of the plot is displyed, by default True
        axes : bool, optional
            if True the axes of the plot are shown, by default True
        colorbar : bool, optional
            if Ture the colorbar is shown, by default True
        unit : str, optional
            the unit of the colobar ('mm' or 'microns'), by default 'mm'
        cbarorientation : str, optional
            The orientation of the colorbar ('vertical' or 'horizontal'), by default 'vertical'
        scalars : np.array, optional
            another array with the same dimension that can be used for assigning a new color coding
            to the mesh.
        save : bool, optional
            if True the resultin plot is saved to disk.

        Returns
        -------
        mlab
            mayavi mlab object

        """
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
                print(X1.shape, Y1.shape, array.shape)
                mesh = mlab.mesh(X1[:array.shape[0],:array.shape[1]],
                                 Y1[:array.shape[0],:array.shape[1]], array,
                                 scalars=scalars, mask=array.mask)
            else:
                print(X1.shape, Y1.shape, array.shape)
                mesh = mlab.mesh(X1[:array.shape[0],:array.shape[1]],
                                 Y1[:array.shape[0],:array.shape[1]],
                             array, scalars=scalars)

        else:
            print(X1.shape, Y1.shape, array.shape)
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
        """3D plot of the distance array using matplotlib.

        This function uses matplotlib for producing a 3D plot of the distances array.
        Matplotlib is slower compared to Mayavi hence is preferable to use self.tredplotmayavimesh.

        Parameters
        ----------
        arr : np.array, optional
            if None self.array is used, by default None
        stride : str, optional
            the strides determine the sampling of the original array strides to 1 corresponds to use all the data, by default 'auto'
        plot : bool, optional
            if Ture the result is shown, by default True
        title : str, optional
            the title of the plot, by default '3D Plot'
        zprop : bool, optional
            proportion of the z axis if True the aspect ratio of the axes is preserved, by default True
        cm : str, optional
            The colormar used for the colorbar, by default 'cw'

        Returns
        -------
        matploltlib.plt
            matplotlib.plt object
        """
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
                print('Change stride to:', stride)
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

    def delcol(self, v=False, col2remove=-1):
        """Deletes columns form the array

        It takes the name of the output file form the profilometer. It returns
        an array with the data rearranged. Deleting the last number of the row,
        this can solve the issue of misalingment of the rows.

        Parameters
        ----------
        v : bool, optional
            verbose mode for diagnostic purposes, by default False
        col2remove : int, optional
            the number of column to be removed 0 for the first column, by default -1
        """
        self.inizialize()
        b = np.delete(self.array, np.s_[col2remove], axis=1)
        rowtoflip = list(range(1, self.parameters.numrows, 2))
        b[rowtoflip, :] = np.fliplr(b[rowtoflip, :])  # sobstitute rows
        self.diagnose(b, self.parameters.numcols - 1,
                      self.parameters.numrows, v)  # Run diagnose errors
        self.array = b
        self.log['method'] = 'delcol (for deleting a column)'
        self.log['removed col'] = col2remove

    def fliprows(self, snr='.dist', v=False):
        # TODO Update this function with the results of the thesis
        '''
        It takes the name of the output file form the profilometer. It returns
        an array with the data rearranged.


        '''
        rowtoflip = list(range(1, self.parameters.numrows, 2))
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
        """Detect the corners of a rectangular sample

        This method uses OpenCV template matching algorithms for finding the corners of the samples.
        The coordinates of the corners are stored in self.corners.
        For performing the template matching the mask array is used by default, the user must acquire the sample
        so that the background will be filtered by the self.mfilter function.
        For templates are created representing the corners of the sample.
        Before using the template matching algorithm a morphological open and closed is performed for improving the results
        of the matching.
        For improving the performance, the original surface can be devide in four region, each region
        containing a corner so that the template matching can be performed only on the region where the
        corner should be found. In this case the sample must be centered at least approximately.

        Parameters
        ----------
        data : np.array, optional
            the array to be used for performing the corner detection by default is used the mask of the
            measurment that shows good corners if the sample is in range and the background is out of range, by default 'mask'
        w : int, optional
            width of the templated to be used for corner detection, by default 100
        h : int, optional
            height of the template to be used for corner detection, by default 100
        Kernel_1 : tuple, optional
            dimension of the kernel for performing a preprocessing of the data, by default (15,15)
        Kernel_2 : tuple, optional
            dimension of the kernerl for the second filter for performina preprocessin of the data, by default (200,200)
        pad_width : int, optional
            the dimension of the pixels added to the scan for allowing the identification of corners very close to the borders, by default 100
        subdivide : bool, optional
            if Ture the sample is subdivided so that the template is computed only on one part of the sample, by default True
        method : str, optional
            The method used for computing the score of the matching of the template, by default "cv2.TM_CCORR_NORMED"
        plot : bool, optional
            if Ture a plot of the results is shown at the end, by default True

        Raises
        ------
        Warning
            Template widht and height should be even.
        """
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
        t_TL[h//2:,w//2:] = 0 #Top right
        t_TR[h//2:,:w//2] = 0 #Top left
        t_BL[:h//2,w//2:] = 0 # Bottom right
        t_BR[:h//2,:w//2] = 0 #Bottom left

        templates = {"TR":t_TR, "TL":t_TL,"BR": t_BR,"BL": t_BL}

        # Image preprocessing
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
            # we subdivide the surface in 4 areas and we search for the template
            # only in the are wher it should be
            if subdivide:
                if j == 'TL':
                    img = img[:hz//2,:wz//2]
                if j == 'TR':
                    wp = wz//2
                    img = img[:hz//2,wz//2:]
                if j == 'BL':
                    hp = hz//2
                    img = img[hz//2:,:wz//2]
                if j == 'BR':
                    wp = wz//2
                    hp = hz//2
                    img = img[hz//2:,wz//2:]
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

            corner = (int(top_left[0] + w//2 -pad_width + wp),
                      int(top_left[1] + h//2 -pad_width + hp))
            self.corners[j] = corner
            corners.append(corner)


        # PLOTTING
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
        '''Correction of the centroids
        This function re-align the centroids and create a new matrix.
        When the probe is acquiring in dynamic mode  it acquires signal
        while is running. If a snake pattern is used mantaining the same  starting
        and ending point for each row or column the real center of the area that has
        been analyzed by the laser is  shifted toward the scanning direction
        of a distance that we call integration distance (see fig. a.).
        This function correct this shift using spline, creating a matrix
        Xc (cfr. fig. a) interpolating the value from the array for every line
        and recopmputing the values for an allingn matrix (cfr. fig. b).

        In many systems this shifting is corrected during the acquistion.

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
        """Subtract a plane from the array of distances.


        See https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6
        for further documentation.

        Parameters
        ----------
        plot : bool, optional
            if Ture the plane subtracted is plotted, by default False
        inverted : bool, optional
            if True the resulting array is inverted (peaks are positive values), by default True
        eliminate_mask : bool, optional
            if True the mask of the original array is deleted, by default False
        dil : int, optional
            the plane plotted is diluted by the factor dil allowing a faster rendering (this does not affect
            ) the computation, by default 10
        array : np.array, optional
            the data to be used for performing the plane subtraction if None self.array is used, by default None
        order : int, optional
            the order of the polynomial fit if 1 a plane is fitted if 2 a surface , by default 1
        zprop : bool, optional
            allows to keep the right aspect ratio on the plot, by default False
        refplane : bool, optional
            if True the reference plane is calculated, by default True
        QC : bool, optional
            if True the function is used for assesing the quality of the measurment, by default False
        returncalculated_plane : bool, optional
            if True the calculated plane is returned, by default False

        Returns
        -------
        np.array
            returns the distance matrix after the subtraction or the reference plane calculated (if
            returncalculated_plane is True.)
        """

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
                z = self.array[~array.mask].data.flatten()

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


        print(r"X axis inclination  %s degrees" %(xin))
        print(r"Y axis inclination  %s degrees" %(yin))

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
                             useLens = False,x=None,y=None):
        """Select manually a set of points to be used for fitting a plane and subtracting it
        from the distance array.

        Parameters
        ----------
        plot : bool, optional
            if True the plane subtracted is plotted, by default False
        inverted : bool, optional
            if True the array is subtracted to the plane (in this way the positive values are peaks), by default True
        dil : int, optional
            the dilution factor for speeding up the rendering of the plot, by default 4
        useLens : bool, optional
            if True the working range of the lenses is used for visualizing the original matrix, by default False
        x : array, optional
            the x coordiantes of the points to be used for subtracting the plane, by default None
        y : array, optional
            the y coordiantes of the points to be used for subtracting the plane, by default None

        Returns
        -------
        np.array
            the coordiantes of the points used for the computation as a numpy array.
        """
        import scipy.linalg
        from mpl_toolkits.mplot3d import Axes3D

        # regular grid covering the domain of the data Xr, Yr are undersampled
        # plane  to be plot faster
        X, Y = np.meshgrid(np.arange(0, self.array.shape[1], 1),
                           np.arange(0, self.array.shape[0], 1))


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

        if x is None and y is None:
            while True:
                pts = np.asarray(plt.ginput(-1, timeout=-1)).astype(int)
                if plt.waitforbuttonpress():
                    break
            x = pts[:,0]
            y = pts[:,1]
        else:
            pts = zip(x,y)
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
        return pts

    def sigmafilter(self, sigmas=1):
        """A filter to mask values grater than a certian number of standard deviation (sigmas)

        Parameters
        ----------
        sigmas : int, optional
            Number of standard deviation to be used for filtering, by default 1
        """
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
        """Mask filter apply a mask to the values that don't meet some quality criteria.

        This is the mask filter that can be used to avoid bad values to be take into
        calulation it affects function as subtract bestplane, and plot every
        numpy function used will not take in account the masked value.

        It can maks also value that have a SNR less then argument snrs (default
        is 500 as suggested by Optimet).

        snr can also take snr matrix as input otherwise it will try to use the
        SNR matrix in the folder.

        Parameters
        ----------
        minz : float, optional
            [description], by default None
        maxz : float, optional
            values greater than this distance will be masked if None the nominal working range of the lens is used, by default None
        badto : str, optional
            invalid value can be changed to a specific values this can be 'zero','max' (maximum value), 'none', or np.nan by default False
        add : int, optional
            add an offset to minz and maxz, by default 0
        snr : bool, optional
            if Ture the snr matrix is used in the computation, by default True
        lens : bool, optional
            if True and minz and maxz are None the working range of the lens is used, by default True
        total : str, optional
            can be used for setting the threshold for the total filtering ('aceptable' or 'optimal'), by default 'optimal'
        snrs : int, optional
            the threshold for the SNR filter, by default 512

        Returns
        -------
        [type]
            [description]
        """
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

        if minz != None and maxz != None:
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
            rowtoflip = list(range(1, self.parameters.numrows, 2))
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
                print("No total file found!")

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
        """Plot the Amplitude Distribution function

        This function shows the distribution of heights as an histogram.
        This is called Amplitude Distribution function. The vertical lines
        represent the working range of the probe if the lens object is available.

        Parameters
        ----------
        vminx : float or str, optional
            the minumim value of the x axis, accept also '-sigma','-2sigma',
            '-3sigma', by default None
        vmaxx : float or str, optional
            the maximum value of the x axis accept also 'sigma','2sigma','3sigma',
             by default None
        auto : bool, optional
            finde automatically the min and max value of the x axis, by default True
        vlines : bool, optional
            if true vertical lines representing the working range if lens object is
            available are shown, by default True
        shift : int, optional
            the shifting of the vertical lines, by default 0
        plot : bool, optional
            if Ture the plot is shown, by default True
        title : bool, optional
            the title of the plot, by default True

        Returns
        -------
        matplotlib.pyplot object
            matpltolib plot object
        """
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
                    linestyles='dashed')
                ax2.vlines(
                    vmaxx + shift,
                    0,
                    ax2.axes.get_ylim()[1],
                    color='r',
                    linestyles='dashed')
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
            print(min(x[1]), max(x[1]))
            print('STDV: %s' % (np.std(np.ravel(self.array))))

            ax2.vlines(
                self.lens.LENS_MINd,
                0,
                ax2.axes.get_ylim()[1],
                color='r',
                linestyles='dashed')
            ax2.vlines(
                self.lens.LENS_MAXd,
                0,
                ax2.axes.get_ylim()[1],
                color='r',
                linestyles='dashed')

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
        """Test and plot the timing diagram of the CCD and the axis trigger

        This function uses the axis velocity and the stage step for calculating
        the trigger frequency. It is mandatory for calculating the exposure lenght.
        Eventually, it is possible to plot the timing diagram, for comparing
        the frequency of the trigger with the frequency of the CCD of the probe.
        Note that if the frequency of the trigger is not a multiple of the
        frequency of the CCD the probe might lose many measuremets, resulting in a
        high number of missing values.

        Parameters
        ----------
        xlim : str, optional
            the limit of the x axis, by default 'auto'
        plot : bool, optional
            if True the plot is shown, by default True

        Returns
        -------
        matplotlib.pyplot
            matplotlib.pyplot object for plotting
        """
        from squarewave import sqwav
        time = (self.parameters.rangeX) / float(self.parameters.X_axis_vel)
        Freq = self.parameters.numcols / time
        Freq2 = 1 / (self.parameters.stage_step /
                     float(self.parameters.X_axis_vel))  # using velocity
        print('Freq. from column:', Freq, 'Frq.:', Freq2, 'CCD:', self.parameters.CCDfreq)
        integrationdistance__mm = (
            1 / self.parameters.CCDfreq) * float(self.parameters.X_axis_vel)
        real_x_resolution = integrationdistance__mm * 1000 + self.lens.X_laser_Spot_Size
        self.log['Calculated spot diameter'] = real_x_resolution
        print('Integration distance: %s mm Signal averaged area: %s micron' % (integrationdistance__mm, real_x_resolution))
        if Freq2 > 3000:
            print("Trigger frequency too high!")
        names = ['CCD', 'Trigger', ]
        # put always the lowest frequency at the end
        freq = [self.parameters.CCDfreq, Freq2, ]
        Triggerduration = 1000 / Freq2 - 0.0005
        durations = ['half period', Triggerduration]
        offsets = [0, 0.0005]
        a = list(zip(freq, durations, offsets, names))
        y_offset = 0
        for i in a:
            sqwav(i[0], i[1], i[2], i[3], y_offset=y_offset, xlim=xlim)
            y_offset += 1.5
        plt.yticks(np.arange(0.5, y_offset, 1.5), names)
        plt.ylim(-0.2, y_offset)
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
        """Plots an estimation of the spot irradiate by the laser and the pattern.

        This function gives you an insight on the resolution, when performing
        acquisition in dynamic mode. It shows you the
        dimension of the laser spot size its initial and final postions, and
        the space that the laser has been moving while recording the line.
        The disalingment due to the CCD frequency can be evaluated.

        '''

        Parameters
        ----------
        centroid : bool, optional
            if True a plot of the centroid is shown by default True
        angle : int, optional
            the angle of the elipse representing the laser spot, by default 0
        laserdirection : str, optional
            The direction of the laser spot (laser diodes might have elipsodial shape), by default 's'
        shifterror : bool, optional
            [description], by default False
        plot : bool, optional
            if True the plot is shown., by default True
        """

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

        XY = list(zip(X2.flatten(), Y2.flatten()))
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

            XY = list(zip(X3.flatten(), Y3.flatten()))
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
        """Select interactively a regions of interest it the acquisition.

        Parameters
        ----------
        scalingfactor : int, optional
            the scalingfactor is used when the whole acquisition can not be plot
            due to memory limitations, by default 1
        useSNR : bool, optional
            instead of showing the distances matrix is uses the SNR matrix, the ROI will be however taken
            from the distance matrix, by default False
        save : str, optional
            if 'ROI' the ROI are save in self.ROI if 'overwrite' in self.array if
            'indexes' only the indexes are saved, by default 'ROI'
        suffix : str, optional
            the suffix to be used when saving the ROI, by default '_ROI'
        """
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
            print("(%s, %s) --> (%s, %s)" % (x1, y1, x2, y2))
            y1=y1*scalingfactor
            y2=y2*scalingfactor
            x1=x1*scalingfactor
            x2=x2* scalingfactor

            if save == 'ROI':
                self.ROI = self.crop(y1=y1,y2=y2,x1=x1,x2=x2,
                                     copy =True,
                                     suffix = suffix)
                print("Last ROI saved in self.ROI")
            if save == 'overwrite':
                self.crop(y1=y1,y2=y2,x1=x1,x2=x2,
                                     copy =False,
                                     )
            if save == 'indexes':
                self.ROIs_indices[suffix+str(counter[0])] = (y1,y2,x1,x2)
                counter[0] +=1
                print("All ROI indices saved in self.ROI_indices")
            else:
                raise ValueError("%s is not a valid parameter for save argument" %save)



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
        """Crop acquisition and recompute the acquisition dimension

        This is a better approach for cropping the array. It allows
        the creation of new ns object with deep copying using the
        copy arguments.

        After the cropping, it recalculates:
        surf.parameters.numcols
        surf.parameters.numrows
        surf.parameters.rangeX
        surf.parameters.rangeY

        Parameters
        ----------
        y1 : int
            the intial y value (corresponding to the row)
        y2 : int
            the final y value (correspoinding to the last row)
        x1 : int
            the intial x value (correspoinding to the column)
        x2 : int
            the final x value (corresponding to the column)
        copy : bool, optional
            if True the cropo does not overwrite self.array, by default False
        suffix : str, optional
            the suffix to be appended to the self.name, by default False

        Returns
        -------
        ns instace
            if copy=True it return a surface object ns mantaining the oringal.
        """

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
        """Crop a region of interest centered on the accquisition.

        This method can be used for cropping the scan to get a square region of interest.
        It finds the center of the acquisition and crop a square with the side equal to
        square_high_mm.
        If usecroners==True it uses the coreners detected using self.cornersdetector.

        Parameters
        ----------
        square_high__mm : int, optional
            side of the square of the ROI, by default 1
        usecorners : bool, optional
            if True the corners found using the self.cornersdetector are used, by default False
        corner : str, optional
            corner to be used for testing (TO BE REMOVED), by default 'TL'
        """
        #calculate the side of square in pts
        pts = int(round(square_high__mm/self.parameters.stage_step/2,0))
        if usecorners:
                if self.corners[corner] != None:
                    x1 = (self.corners['TL'][0] + self.corners['BL'][0])//2
                    x2 = (self.corners['TR'][0] + self.corners['BR'][0])//2
                    y1 = (self.corners['TL'][1] + self.corners['TR'][1])//2
                    y2 = (self.corners['BL'][1] + self.corners['BR'][1])//2

                    xcent = int(round((x2 - x1)/2 + x1,0))
                    ycent = int(round((y2 - y1)/2 + y1,0))

                else:
                    print("No corners found!")
                    return

        else:
            ycent, xcent = np.array(self.array.shape)//2

        self.array = self.array[ycent-pts:ycent+pts,
                                xcent-pts:xcent+pts]
        self.parameters.numcols = y2-y1
        self.parameters.numrows = x2-x1
        self.parameters.rangeX = square_high__mm
        self.parameters.rangeY = square_high__mm



    def concentric_ROI(self,step=20,plot='show',save=False):
        """Subdivide the area in concentric ROI with increasing side
         and compute metrological parameters for each iteration increases
         the side of the ROI by 2*step*i.

         This procedure is helpful for understanding the homogeneity of the ROI.

        Parameters
        ----------
        step : int, optional
            the number of  measurments of HALF of the side of the ROI.
        plot : str, optional
            if 'show' the resulting plot is show, by default 'show'
        save : bool, optional
            if True the results are saved, by default False

        Returns
        -------
        tuple
            a tuple containing four lists with the distancefromcenter,Sq,
            Rskw, Rkurt
        """
        from scipy.stats.mstats import skew, kurtosis
        ycent, xcent = np.array(self.array.shape)//2
        distancefromcenter = []
        Sq=[]
        Rskw=[]
        Rkurt=[]
        rectangles = []

        i=1
        while step*i<xcent and step*i<ycent:
            ROI = self.array[ycent-step*i:ycent+step*i,xcent-step*i:xcent+step*i].flatten()
            print(ROI.size)
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
                ACD.savefig(os.path.join(getcwd(),'Results','%s ConcentricROIPlot.png'%self.name))
                plt.close()



        return distancefromcenter,Sq,Rskw, Rkurt

    def margin_ROI(self,step=20,plot='show',save=False):
        """Compute metrological paramters excluding an increasing portion of the margin

        This method is analogues to self.concentric_ROI but it does not require the sample to be
        squared. For each iteration excludes a greater boarder from the commputation starting from
        the a margin excluding the value assigne to the variable *step* and increasing the margin of
        step*i where i is the iteration number.

        Parameters
        ----------
        step : int, optional
            initial with of the margin excluded by the compuatation in measurments, by default 20
        plot : str, optional
            if 'show' a plot showin the procedure is shown at the end, by default 'show'
        save : bool, optional
            if True the results are save to file, by default False

        Returns
        -------
        tuple
            a tuple containing four lists with the value computed for each iteration
            margin, Sq, Rskw, Rkurt.
        """
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
            print(ROI.size)
            margin.append(step*i)
            Rskw.append(skew(ROI))
            Rkurt.append(kurtosis(ROI))
            Sq.append(np.std(ROI)*1000)
            rectangles.append((step*i,step*i,self.array.shape[1]-step*i*2,
                               self.array.shape[0]-step*i*2))
            print((step*i,step*i,self.array.shape[1]-step*i*2,
                               self.array.shape[0]-step*i*2))
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
                ACD.savefig(os.path.join(getcwd(),'Results','%s MarginsROIPlot.png'%self.name))
                plt.close()



        return margin,Sq,Rskw, Rkurt

    def calculate_profile_metrology(self,profile=None,roundx=3):
        """Calculates some amplituted metrology parameters for a profile.

        This function was meant to be used for computing the surface metrolgoy
        parameters for the profiles using a discrete aproximation.

        # TODO: define sampling lenght and continuous calculation

        Parameters
        ----------
        profile : np.array, optional
            The profile to be used, by default None
        roundx : int, optional
            The number of digit to be rounded, by default 3
        Returns
        -------
        dict
            a dictionary structured as followed
            {'Rz':Rz, 'Rv':Rv, 'Rt':Rt, 'Ra':Ra, 'Rp':Rp, 'Rq':Rq,
                        'Rsk':Rsk, 'Rku':Rku, 'Rc':Rc}
        """
        ##  DEFINE SAMPLING LENGHT !!!!
    ####
        from scipy.stats.mstats import skew, kurtosis
        if self.log['Best plane correction']== 'Not Performed':
            print('WARINING BEST PLANE NOT PERFORMED!')
        if profile == None:
            data=self.array[0]
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
        """Calculates area metrology parameters.

        Parameters
        ----------
        ROI : [type], optional
            uses the indexes of in the list for selecting a region of interest
            data=self.array[ROI[0]:ROI[1],ROI[2]:ROI[3]], by default None
        roundx : int, optional
            [description], by default 3
        margin : int, optional
            a margin can be applied for excluding measurments in the boarder, by default None
        subtractplane : bool, optional
            if True a plane subtraction is performed on the data (useful if we are using a ROI)
            and we did not subtract the plane before, by default False

        Returns
        -------
        dict
            a dictionary structured as followed
            dicmetrology = {'Sz_microns':Sz, 'Sv_microns':Sv, 'St_microns':St,
            'Sar_microns':Sar, 'Sp_microns':Sp, 'Sq_microns':Sq,
                        'Ssk':Ssk, 'Sku':Sku, 'Sa_microns':Sa}
        """
        from scipy.stats.mstats import skew, kurtosis
        if self.log['Best plane correction'] == 'Not Performed':
            print('WARINING BEST PLANE NOT PERFORMED!')
        else:
            if subtractplane and self.log['Best plane correction']== 'Performed':
                print('PLANE ALREADY SUBTRACTED!')
            else:
                print('calculating...')

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

    def go(self):
        """Perform self.mfilter self.subtractplane and self.plot with the default parameters.
        """
        self.mfilter()
        self.subtractplane()
        self.plot()

# Exporting tools.
    def GEOTIFF(self):
        """ Export the data in GEOTIFF using gdal

        Save a raster from a C order array.

        +proj=tmerc +lat_0=51.4 +lon_0=7 +k=1 +x_0=0 +y_0=0 +ellps=WGS84
        +towgs84=0,0,0,0,0,0,0 +units=m +no_defs

        """
        import osgeo.gdal as gdal
        self._mk_results_folder()

        dst_filename = path.join("Results","%s.tiff" %self.name)
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
        """Export of the flatten self.array parameter.
        """
        self._mk_results_folder()

        with open(path.join('Results', '%s.sdist' %self.name), "wb") as f:
            newFileByteArray = bytearray(self.array.flatten())
            f.write(newFileByteArray)

    def npy_export(self):
        """Export in .npy format, can be load with np.load()
        """
        np.save(path.join('Results',self.name), self.array)
        self.savelog()

    def mat_export(self, mdict = None, savemask =False):
        """Export in .mat format.

        Parameters
        ----------
        mdict : dict, optional
            the dictionary to be save can be created manually when
            None is created automatically, by default None
        savemask : bool, optional
            if Ture the mask is save as a separate matrix in the mat dict, by default False
        """
        import scipy.io

        self._mk_results_folder()
        if mdict == None:
            if savemask:
                mdict = {'data' + self.name: self.array.data,
                         'mask' + self.name: self.array.mask}
            else:
                mdict = {'m' +
                    self.name: self.array}
        scipy.io.savemat(path.join('Results',self.name+'.mat'),
            mdict=mdict)
        print("Saved %s in %s" %(self.name+'.mat',path.join(getcwd(),'Results')))
        self.savelog()

    def h5f_export(self, dataset=''):
        """Export in .hdf5 format

        Parameters
        ----------
        dataset : str, optional
            the title of the dataset when '' the self.name is used, by default ''
        """
        import h5py
        if dataset == '':
            dataset = self.name
        h5f = h5py.File(os.path('Results',dataset + '.hdf'), 'w')
        h5f.create_dataset(dataset, data=self.array)
        self.savelog()

    def gwy_export(self,unit = "mm"):
        """Export the data using gwyfile module for Gwiddion

        Parameters
        ----------
        unit : str, optional
            the unit of measurment, by default "mm"
        """
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

        with open(path.join(r'Results', self.name + 'v2' + '.spm'), 'w') as o:
            for i in header:
                o.write(i)

    def SPM_export(self):
        """export the data in SPM format
        """
        self._mk_results_folder()

        n = '\n'
        with open(path.join('Results',self.name + '.spm'), 'w') as o:
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

        print('X1 shape: ', X1.shape, 'Y1 shape: ', Y1.shape)
        X = X1.flatten()
        print('X1:', X.size)
        Y = Y1.flatten()
        print('y1:', Y.size)
        Z = self.array.flatten()
        print('z:', Z.size)
        if invertvalue:
            Z = (self.lens.LENS_MAXd + 3) - Z
        xyz = list(zip(X, Y, Z))
        number_of_vertex = len(Z)
        return xyz, number_of_vertex

    def ASCII_export(self, kind='ply', invertvalue=False, centroids=False):
        """
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
        """

        n = '\n'
        xyz,nv = self.pointCloud_generator(centroids=centroids,invertvalue=invertvalue)
        self._mk_results_folder()
        if kind == 'ply':
            with open(path.join(getcwd(),r'Results',self.name + '.ply'), 'w') as o:
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
            with open(path.join(getcwd(),r'Results',self.name + '.XYZ'), 'w') as o:
                o.write('# Channel: Gray' + n)
                o.write('# Lateral units: m' + n)
                o.write('# Value units: m' + n)
                for i in xyz:
                    o.write('%.4e\t%.4e\t%.4e' % (i[0], i[1], i[2]) + n)

    def exp_appyreport(self, tpath, metrology=3,flipyax=None):
        """Export a report in odt if a tempalte is available

        Parameters
        ----------
        tpath: str,
            the path of the template
        metrology : int, optional
            the number of ROI to be selected in the report, by default 3
        flipyax : int, optional
            if -1 the axis is flipped, by default None
        """
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

        print('Plotting...')
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
            print('Size exceed memory limits,diluting data...')
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
        print('ADF...')
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
        print('Attributes...')
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
            print('Sample image..')
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
        print('SNR...')
        self.SNR(plot=True)[0]
        figfile2 = BytesIO()
        try:
            if flipyax == -1:
                plt.gca().invert_yaxis()
                plt.gca().invert_xaxis()
            plt.savefig(figfile2, format='png', dpi=50)

        except MemoryError:
            print('Size exceed memory limits,diluting data...')
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
        renderer = Renderer(os.join(tpath,"Template_new_metro%s.odt" %(metrology),
                            locals(), '%s.odt' % (self.path)))
        renderer.run()
        #restore original figsize
        plt.rcParams['figure.figsize'] =  originalfigsize
        gc.collect()

    def viewer(self, absolute=False,complete=True):
        """A viewr based on PythonQwt

        Based on https://pypi.python.org/pypi/guiqwt this
        viewer allows to inspect the data.

        Parameters
        ----------
        absolute : bool, optional
            if True absolute coordante systems is used, by default False
        complete : bool, optional
            if Ture also SNR and Total are loaded, by default True

        Returns
        -------
        win
            a PythonQwt object
        """
        # Based on PythonQwt

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
                rowtoflip = list(range(1, self.parameters.numrows, 2))
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
            print("masked data")
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
        """Normalize the data to the range 0-255

        Parameters
        ----------
        astype : type, optional
            the type of the new data, by default np.uint8

        Returns
        -------
        np.array
            The normalized np.array
        """
        return (255*(self.array - np.min(self.array))/np.ptp(self.array)).astype(astype)

    def Segmentation(self, side, kind=None, mfilter=True,
                      overwrite=False,margin=None,norm=False):
        """Divides the scan in squeres with equal side and compute metrological parameters for each squere.

        Parameters
        ----------
        side : int
            the side of the squares.
        kind : int, optional
            the process to be carried out for each square 1 - Roughenss
            2 - Skwenness, 3 - Kurtosis,4 - False Color Composite,
            5 - Equalized False Color Composite, if None the user is aske to chose, by default None
        mfilter : bool, optional
            if True the mfilter is used, by default True
        overwrite : bool, optional
            if True the self.array is overwrited with the results of self.Segmentation, by default False
        margin : int optional
            the margins to be exclueded from the computation in number of measurments, by default None
        norm : bool, optional
            if True the results is divided by the mean if 'to1' is also normalized to 1, by default False

        Returns
        -------
        np.array
            the array with the computed values
        """
        from scipy.stats import skew, kurtosis
        if kind is None:
            print(" Select operation to perform on patches:")
            print(" 1 - Roughenss")
            print(" 2 - Skwenness")
            print(" 3 - Kurtosis")
            print(" 4 - False Color Composite")
            print(" 5 - Equalized False Color Composite")
            kind = eval(input('Type the number: '))


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
            print('Eliminated %i columns and %i row' % (restoc, restor))
        elif restoc != 0:
            final = filteredarr[::, :-restoc]
            print('Eliminated %i columns ' % (restoc))
        elif restor != 0:
            final = filteredarr[:-restor, ::]
            print('Eliminated  %i row' % (restor))

        else:
            final = filteredarr

        numcolsa = int(final.shape[1] / side)
        print(numcolsa)
        numrowsa = int(final.shape[0] / side)
        print(numrowsa)

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
                path.join('Results','RGB composite %s pixel%s.tiff' %
                (self.name, side)))
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
                print('Divinding by mean...')
                if norm =='to1':
                    minval = np.nanmin(self.array)
                    print('Normalizing...')
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
        """Divides the scan area in contiguous squares of increasing area
        and compute metrological parameters for each step

        This function iterates self.Segmentation for square of different side
        the computation carried out for each square is the same of
        self.Segmentation.
        1 - Std. Deviation Roughness
        2 - Skewness
        3 - Kurtosis
        6 - This is a 2D implementation of the 1D algorithm proposed in:
        The multi-scale roughness analyses and modeling of abrasion with the grit
        size effect on ground surfaces.
        Wear 286-287 (2012) 124-135  Maxence Bigerelle et al.

        Parameters
        ----------
        kind : int
            the kind of computation required
        step : int
            the increasing factor of the square side
        plot : bool, optional
            if True the results are plotted, by default True
        save : bool, optional
            if True the results are saved to disk, by default True
        xlim : float, optional
            the xlim of the plot, by default None
        ylim : tuple, optional
            the ylim of the plot, by default None
        norm : bool, optional
            if Ture the results are divided by the mean if 'to1' are also normalized to 1, by default False

        Returns
        -------
        dict
            a dictionary with the results of self.Segmentation for each step
        """
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
            steps = list(range(2,maxside,step))
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
        """
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

        """
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

        calcRoughnessANS = input('Do you want to calulate roughness? y/n ')
        if calcRoughnessANS.lower() == 'y':
            for i in np.unique(segments_slic):
                segment = np.ma.masked_where(segments_slic != i, self.array)
                Roughness = np.std(segment)
                self.array[segment.nonzero()] = Roughness

        self.log['Processing'] = self.log['Processing'].join(
            'SLIC minz=%s,maxz=%s,n_segments=%s, compactness=%s, enforce_connectivity=%s' %
            (minz, maxz, n_segments, compactness, enforce_connectivity))

        calccraksANS = input('Do you want to calulate craks? y/n ')

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
        """This is an utility function that creates the folder results if does not exists
        """
        if not path.exists('Results'):
            makedirs(path.join(getcwd(),'Results'))

    def _save_plot(self,fig,name,dpi=300):
        """Save the plot to disk

        This is an utility function for saving the
        plot to disk.

        Parameters
        ----------
        fig : fig
            a matplotlib figure type object
        name : str
            the name of the figure
        dpi : int, optional
            the dot per inch of the figure
        """
        self._mk_results_folder()
        while True:
            try:
                fig.savefig(path.join('Results',name), dpi = dpi)
                break
            except:
                dpi /= 2
        plt.close(fig)
