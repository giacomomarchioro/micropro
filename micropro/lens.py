'''
This class contains information regarding the lens of the conoprobe.
'''
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