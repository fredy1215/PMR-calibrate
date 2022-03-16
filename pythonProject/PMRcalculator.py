from PyQt5 import QtCore, QtGui, QtWidgets
import scipy.io as sco
import os
import fnmatch
import numpy as np
from scipy.optimize import fsolve


def find_files(filename, search_path):
    result = []

    if filename.find('*') != -1:

        for root, dir, files in os.walk(search_path):
            for a in files:
                if fnmatch.fnmatch(a, filename):
                    result.append(os.path.join(root, a))
        return result
    else:
        for root, dir, files in os.walk(search_path):
            if filename in files:
                result.append(os.path.join(root, filename))
        return result


File = find_files('*.mat', 'C:\\Users\\Meng Lu\\Box\\Thesis\\code\\app')
matstruct_contents_1_5 = sco.loadmat(File[0])
matstruct_contents_3 = sco.loadmat(File[1])


class Ui_PMRcalculater(object):
    def RDselect(self):
        if self.RB1_5.isChecked() == 1:
            self.SPT1myoT.setText('995.8')
            self.TI = matstruct_contents_1_5['TI']
        elif self.RB3.isChecked() == 1:
            self.SPT1myoT.setText('1229')
            self.TI = matstruct_contents_3['TI']

    def TIpushed(self):
        H = list(range(40, 100 + 1, 1))
        Nflash = list(range(20, 50 + 1, 1))
        FA = list(range(10, 20 + 1, 1))
        Echo = list(np.arange(3, 6 + 0.1, 0.1))
        Echo = [round(elem, 4) for elem in Echo]
        if self.HrT.text() == '':
            msgTIpushed = QtWidgets.QMessageBox()
            msgTIpushed.setWindowTitle('error(TI)')
            msgTIpushed.setText('Please enter all the fields in step 1!')
            msgTIpushed.exec_()
        else:
            A = H.index(int(self.HrT.text()))
            B = Nflash.index(int(self.SegT.text()))
            C = FA.index(int(self.FaT.text()))
            D = Echo.index(float(self.EchosT.text()))
            TIselected = self.TI[A][B][C][D]
            self.TiT.setText(str(TIselected))

    def CalT1pushed(self):
        if self.ScannerPMRT.text() == '':
            msgT1pushed = QtWidgets.QMessageBox()
            msgT1pushed.setWindowTitle('error(T1)')
            msgT1pushed.setText('Please enter all the fields in step 1 and 2!')
            msgT1pushed.exec_()
        else:
            TI = float(self.TiT.text())
            echos = float(self.EchosT.text())
            n = float(self.SegT.text())
            HR = float(self.HrT.text())
            TR = 1 / HR * 60 * 1000
            flipA = float(self.FaT.text())
            M0 = 1
            PMR = float(self.ScannerPMRT.text())
            T1_myo = float(self.SPT1myoT.text())

            def eq(T1_pl):
                return (abs(M0 * (np.exp(-TI / T1_pl) - 1) - (np.exp(-TI / T1_pl) * (M0 * (np.exp((echos * n - TR + TI)
                                                                                                  / T1_pl) - 1) + np.exp(
                    (echos * n - TR + TI) / T1_pl) * (M0 * np.exp(-(echos * n) / T1_pl) * np.cos((np.pi * flipA)
                                                                                                 / 180) ** n * (
                                                              np.exp(-TI / T1_pl) - 1) + (
                                                              M0 * (np.exp(-echos / T1_pl) - 1) * ((np.exp(
                                                          -echos / T1_pl) * np.cos((np.pi * flipA) / 180)) ** n - 1))
                                                      / (np.exp(-echos / T1_pl) * np.cos(
                            (np.pi * flipA) / 180) - 1)))) / (
                                    np.exp(-(echos * n) / T1_pl) * np.exp(-TI / T1_pl) * np.exp((echos * n - TR + TI)
                                                                                                / T1_pl) * np.cos(
                                (np.pi * flipA) / 180) ** n + 1)) / abs(
                    M0 * (np.exp(-TI / T1_myo) - 1) - (np.exp(-TI / T1_myo) * (M0 * (np.exp((echos * n - TR + TI)
                                                                                            / T1_myo) - 1) + np.exp(
                        (echos * n - TR + TI) / T1_myo) * (M0 * np.exp(-(echos * n) / T1_myo) * np.cos(
                        (np.pi * flipA) / 180) ** n * (np.exp(-TI / T1_myo)
                                                       - 1) + (M0 * (np.exp(-echos / T1_myo) - 1) * (
                            (np.exp(-echos / T1_myo) * np.cos((np.pi * flipA) / 180)) ** n - 1)) / (
                                                                   np.exp(-echos / T1_myo)
                                                                   * np.cos((np.pi * flipA) / 180) - 1)))) / (
                            np.exp(-(echos * n) / T1_myo) * np.exp(-TI / T1_myo) * np.exp(
                        (echos * n - TR + TI) / T1_myo)
                            * np.cos((np.pi * flipA) / 180) ** n + 1)) - PMR)

            # T1_pl = fsolve(eq,[1000,500,1500])
            T1_pl = fsolve(eq, [200, 400, 600, 800, 1000, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])
            T1_plfloor = np.floor(T1_pl)
            max_count = (0,0)
            for num in T1_plfloor:
                appearance = np.sum(T1_plfloor==num)
                if appearance > max_count[0]:
                    max_count = (appearance,num)
            self.PT1T.setText('%.4f' % max_count[1])
                #if T1_pl[i] > 0 and T1_pl[i] <= 3000:
                 #   self.PT1T.setText('%.4f' % T1_pl[i])
                  #  break

            # self.PT1T.setValidator(self.onlyfloat)
            # if float(self.ScannerPMRT.text())>2.5:
            # self.PT1T.setText('%.4f' % T1_pl[1])

            # elif float(self.ScannerPMRT.text())<1:
            # self.PT1T.setText('%.4f' % T1_pl[2])
            # else:
            # self.PT1T.setText('%.4f' % T1_pl[0])

    def CalPMRpushed(self):
        if self.PT1T.text() == '':
            msgTIpushed = QtWidgets.QMessageBox()
            msgTIpushed.setWindowTitle('error(PMR)')
            msgTIpushed.setText('Please enter all the fields in step 1 2 and 3!')
            msgTIpushed.exec_()
        else:
            TI = float(self.SPInversionTimeT.text())
            echos = float(self.SPEchosT.text())
            n = float(self.SPSegT.text())
            HR = float(self.SPHrT.text())
            TR = 1 / HR * 60 * 1000
            flipA = float(self.SPFaT.text())
            M0 = 1;
            T1_pl = float(self.PT1T.text())
            T1_myo = float(self.SPT1myoT.text())

            # def eq(PMR):
            # return (abs(M0 * (np.exp(-TI / T1_pl) - 1) - (np.exp(-TI / T1_pl) * (M0 * (np.exp((echos * n - TR + TI)
            #                                                                                                  / T1_pl) - 1) + np.exp(
            #                    (echos * n - TR + TI) / T1_pl) * (M0 * np.exp(-(echos * n) / T1_pl) * np.cos((np.pi * flipA)
            #                                                                                                / 180) ** n * (
            #                                                              np.exp(-TI / T1_pl) - 1) + (
            #                                                              M0 * (np.exp(-echos / T1_pl) - 1) * ((np.exp(
            #                                                          -echos / T1_pl) * np.cos((np.pi * flipA) / 180)) ** n - 1))
            #                                                     / (np.exp(-echos / T1_pl) * np.cos(
            #                           (np.pi * flipA) / 180) - 1)))) / (
            #                                   np.exp(-(echos * n) / T1_pl) * np.exp(-TI / T1_pl) * np.exp((echos * n - TR + TI)
            #                                                                                               / T1_pl) * np.cos(
            #                               (np.pi * flipA) / 180) ** n + 1)) / abs(
            #                  M0 * (np.exp(-TI / T1_myo) - 1) - (np.exp(-TI / T1_myo) * (M0 * (np.exp((echos * n - TR + TI)
            #                                                                                           / T1_myo) - 1) + np.exp(
            #                       (echos * n - TR + TI) / T1_myo) * (M0 * np.exp(-(echos * n) / T1_myo) * np.cos(
            #                       (np.pi * flipA) / 180) ** n * (np.exp(-TI / T1_myo)
            #                                                     - 1) + (M0 * (np.exp(-echos / T1_myo) - 1) * (
            #                          (np.exp(-echos / T1_myo) * np.cos((np.pi * flipA) / 180)) ** n - 1)) / (
            #                                                                 np.exp(-echos / T1_myo)
            #                                                                 * np.cos((np.pi * flipA) / 180) - 1)))) / (
            #                          np.exp(-(echos * n) / T1_myo) * np.exp(-TI / T1_myo) * np.exp(
            #                       (echos * n - TR + TI) / T1_myo)
            #                           * np.cos((np.pi * flipA) / 180) ** n + 1)) - PMR)
            PMR = abs(M0 * (np.exp(-TI / T1_pl) - 1) - (np.exp(-TI / T1_pl) * (M0 * (np.exp((echos * n - TR + TI)
                                                                                            / T1_pl) - 1) + np.exp(
                (echos * n - TR + TI) / T1_pl) * (M0 * np.exp(-(echos * n) / T1_pl) * np.cos((np.pi * flipA)
                                                                                             / 180) ** n * (
                                                          np.exp(-TI / T1_pl) - 1) + (
                                                          M0 * (np.exp(-echos / T1_pl) - 1) * ((np.exp(
                                                      -echos / T1_pl) * np.cos((np.pi * flipA) / 180)) ** n - 1))
                                                  / (np.exp(-echos / T1_pl) * np.cos(
                        (np.pi * flipA) / 180) - 1)))) / (
                              np.exp(-(echos * n) / T1_pl) * np.exp(-TI / T1_pl) * np.exp((echos * n - TR + TI)
                                                                                          / T1_pl) * np.cos(
                          (np.pi * flipA) / 180) ** n + 1)) / abs(
                M0 * (np.exp(-TI / T1_myo) - 1) - (np.exp(-TI / T1_myo) * (M0 * (np.exp((echos * n - TR + TI)
                                                                                        / T1_myo) - 1) + np.exp(
                    (echos * n - TR + TI) / T1_myo) * (M0 * np.exp(-(echos * n) / T1_myo) * np.cos(
                    (np.pi * flipA) / 180) ** n * (np.exp(-TI / T1_myo)
                                                   - 1) + (M0 * (np.exp(-echos / T1_myo) - 1) * (
                        (np.exp(-echos / T1_myo) * np.cos((np.pi * flipA) / 180)) ** n - 1)) / (
                                                               np.exp(-echos / T1_myo)
                                                               * np.cos((np.pi * flipA) / 180) - 1)))) / (
                        np.exp(-(echos * n) / T1_myo) * np.exp(-TI / T1_myo) * np.exp(
                    (echos * n - TR + TI) / T1_myo)
                        * np.cos((np.pi * flipA) / 180) ** n + 1))
            self.StdPMRT.setValidator(self.onlyfloat)
            self.StdPMRT.setText('%.4f' % PMR)

    def setupUi(self, PMRcalculater):
        PMRcalculater.setObjectName("PMRcalculater")
        PMRcalculater.resize(778, 496)
        PMRcalculater.setUnifiedTitleAndToolBarOnMac(False)
        self.onlyfloat = QtGui.QDoubleValidator()
        self.onlyfloat.setDecimals(4)
        self.widget = QtWidgets.QWidget(PMRcalculater)
        self.widget.setObjectName("widget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.widget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.gridLayout_7 = QtWidgets.QGridLayout()
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.gridLayout_5 = QtWidgets.QGridLayout()
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.groupBox = QtWidgets.QGroupBox(self.widget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.splitter = QtWidgets.QSplitter(self.groupBox)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName("splitter")
        self.RB1_5 = QtWidgets.QRadioButton(self.splitter)
        self.RB1_5.setObjectName("RB1_5")
        self.RB3 = QtWidgets.QRadioButton(self.splitter)
        self.RB3.setObjectName("RB3")
        self.gridLayout.addWidget(self.splitter, 0, 0, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox, 0, 0, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.SPEchosT = QtWidgets.QLineEdit(self.groupBox_2)
        self.SPEchosT.setObjectName("SPEchosT")
        self.SPEchosT.setValidator(self.onlyfloat)
        self.SPEchosT.setText('5')
        self.gridLayout_2.addWidget(self.SPEchosT, 2, 1, 1, 1)
        self.SPInversionTimeL = QtWidgets.QLabel(self.groupBox_2)
        self.SPInversionTimeL.setObjectName("SPInversionTimeL")
        self.gridLayout_2.addWidget(self.SPInversionTimeL, 0, 0, 1, 1)
        self.SPSegL = QtWidgets.QLabel(self.groupBox_2)
        self.SPSegL.setObjectName("SPSegL")
        self.gridLayout_2.addWidget(self.SPSegL, 1, 0, 1, 1)
        self.SPEchosL = QtWidgets.QLabel(self.groupBox_2)
        self.SPEchosL.setObjectName("SPEchosL")
        self.gridLayout_2.addWidget(self.SPEchosL, 2, 0, 1, 1)
        self.SPFaL = QtWidgets.QLabel(self.groupBox_2)
        self.SPFaL.setObjectName("SPFaL")
        self.gridLayout_2.addWidget(self.SPFaL, 1, 2, 1, 1)
        self.SPInversionTimeT = QtWidgets.QLineEdit(self.groupBox_2)
        self.SPInversionTimeT.setObjectName("SPInversionTimeT")
        self.SPInversionTimeT.setValidator(self.onlyfloat)
        self.SPInversionTimeT.setText('650')
        self.gridLayout_2.addWidget(self.SPInversionTimeT, 0, 1, 1, 1)
        self.SPFaT = QtWidgets.QLineEdit(self.groupBox_2)
        self.SPFaT.setObjectName("SPFaT")
        self.SPFaT.setValidator(self.onlyfloat)
        self.SPFaT.setText('15')
        self.gridLayout_2.addWidget(self.SPFaT, 1, 3, 1, 1)
        self.SPHrT = QtWidgets.QLineEdit(self.groupBox_2)
        self.SPHrT.setObjectName("SPHrT")
        self.SPHrT.setValidator(self.onlyfloat)
        self.SPHrT.setText('65')
        self.gridLayout_2.addWidget(self.SPHrT, 0, 3, 1, 1)
        self.SPT1myoT = QtWidgets.QLineEdit(self.groupBox_2)
        self.SPT1myoT.setObjectName("SPT1myoT")
        self.SPT1myoT.setValidator(self.onlyfloat)
        self.SPT1myoT.setText('995.8')
        self.TI = matstruct_contents_1_5['TI']
        self.gridLayout_2.addWidget(self.SPT1myoT, 2, 3, 1, 1)
        self.SPT1myoL = QtWidgets.QLabel(self.groupBox_2)
        self.SPT1myoL.setObjectName("SPT1myoL")
        self.gridLayout_2.addWidget(self.SPT1myoL, 2, 2, 1, 1)
        self.SPSegT = QtWidgets.QLineEdit(self.groupBox_2)
        self.SPSegT.setObjectName("SPSegT")
        self.SPSegT.setValidator(self.onlyfloat)
        self.SPSegT.setText('30')
        self.gridLayout_2.addWidget(self.SPSegT, 1, 1, 1, 1)
        self.SPHrL = QtWidgets.QLabel(self.groupBox_2)
        self.SPHrL.setObjectName("SPHrL")
        self.gridLayout_2.addWidget(self.SPHrL, 0, 2, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_2, 0, 1, 1, 1)
        self.gridLayout_7.addLayout(self.gridLayout_5, 0, 0, 1, 1)
        self.gridLayout_6 = QtWidgets.QGridLayout()
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.groupBox_3 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_3.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.groupBox_3.setAlignment(QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.groupBox_3.setObjectName("groupBox_3")
        self.formLayout = QtWidgets.QFormLayout(self.groupBox_3)
        self.formLayout.setObjectName("formLayout")
        self.HrL = QtWidgets.QLabel(self.groupBox_3)
        self.HrL.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.HrL.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.HrL.setObjectName("HrL")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.HrL)
        self.HrT = QtWidgets.QLineEdit(self.groupBox_3)
        self.HrT.setObjectName("HrT")
        self.HrT.setValidator(self.onlyfloat)
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.HrT)
        self.SegL = QtWidgets.QLabel(self.groupBox_3)
        self.SegL.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.SegL.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.SegL.setObjectName("SegL")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.SegL)
        self.SegT = QtWidgets.QLineEdit(self.groupBox_3)
        self.SegT.setObjectName("SegT")
        self.SegT.setValidator(self.onlyfloat)
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.SegT)
        self.FaL = QtWidgets.QLabel(self.groupBox_3)
        self.FaL.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.FaL.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.FaL.setObjectName("FaL")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.FaL)
        self.FaT = QtWidgets.QLineEdit(self.groupBox_3)
        self.FaT.setObjectName("FaT")
        self.FaT.setValidator(self.onlyfloat)
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.FaT)
        self.EchosL = QtWidgets.QLabel(self.groupBox_3)
        self.EchosL.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.EchosL.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.EchosL.setObjectName("EchosL")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.EchosL)
        self.EchosT = QtWidgets.QLineEdit(self.groupBox_3)
        self.EchosT.setObjectName("EchosT")
        self.EchosT.setValidator(self.onlyfloat)
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.EchosT)
        self.line = QtWidgets.QFrame(self.groupBox_3)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.SpanningRole, self.line)
        self.TiL = QtWidgets.QLabel(self.groupBox_3)
        self.TiL.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.TiL.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.TiL.setObjectName("TiL")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.TiL)
        self.TiT = QtWidgets.QLineEdit(self.groupBox_3)
        self.TiT.setObjectName("TiT")
        self.TiT.setValidator(self.onlyfloat)
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.TiT)
        self.SuggestTipushButton = QtWidgets.QPushButton(self.groupBox_3)
        self.SuggestTipushButton.setObjectName("SuggestTipushButton")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.SuggestTipushButton)
        self.gridLayout_6.addWidget(self.groupBox_3, 0, 0, 3, 1)
        self.groupBox_4 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_4.setObjectName("groupBox_4")
        self.formLayout_3 = QtWidgets.QFormLayout(self.groupBox_4)
        self.formLayout_3.setObjectName("formLayout_3")
        self.ScannerPMRL = QtWidgets.QLabel(self.groupBox_4)
        self.ScannerPMRL.setObjectName("ScannerPMRL")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.ScannerPMRL)
        self.ScannerPMRT = QtWidgets.QLineEdit(self.groupBox_4)
        self.ScannerPMRT.setObjectName("ScannerPMRT")
        self.ScannerPMRT.setValidator(self.onlyfloat)
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.ScannerPMRT)
        self.gridLayout_6.addWidget(self.groupBox_4, 0, 1, 1, 1)
        self.groupBox_5 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_5.setObjectName("groupBox_5")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_5)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.PT1L = QtWidgets.QLabel(self.groupBox_5)
        self.PT1L.setObjectName("PT1L")
        self.gridLayout_4.addWidget(self.PT1L, 1, 0, 1, 1)
        self.PT1T = QtWidgets.QLineEdit(self.groupBox_5)
        self.PT1T.setObjectName("PT1T")
        self.PT1T.setValidator(self.onlyfloat)
        self.gridLayout_4.addWidget(self.PT1T, 1, 1, 1, 1)
        self.CalT1pushButton = QtWidgets.QPushButton(self.groupBox_5)
        self.CalT1pushButton.setObjectName("CalT1pushButton")
        self.gridLayout_4.addWidget(self.CalT1pushButton, 0, 0, 1, 2)
        self.gridLayout_6.addWidget(self.groupBox_5, 1, 1, 1, 1)
        self.groupBox_6 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_6.setObjectName("groupBox_6")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_6)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.CalStdPMRpushButton = QtWidgets.QPushButton(self.groupBox_6)
        self.CalStdPMRpushButton.setObjectName("CalStdPMRpushButton")
        self.gridLayout_3.addWidget(self.CalStdPMRpushButton, 0, 0, 1, 2)
        self.StdPMRL = QtWidgets.QLabel(self.groupBox_6)
        self.StdPMRL.setObjectName("StdPMRL")
        self.gridLayout_3.addWidget(self.StdPMRL, 1, 0, 1, 1)
        self.StdPMRT = QtWidgets.QLineEdit(self.groupBox_6)
        self.StdPMRT.setObjectName("StdPMRT")
        self.StdPMRT.setValidator(self.onlyfloat)
        self.gridLayout_3.addWidget(self.StdPMRT, 1, 1, 1, 1)
        self.gridLayout_6.addWidget(self.groupBox_6, 2, 1, 1, 1)
        self.gridLayout_7.addLayout(self.gridLayout_6, 1, 0, 1, 1)
        self.horizontalLayout.addLayout(self.gridLayout_7)
        # events
        self.RB1_5.setChecked(1)
        self.RB1_5.toggled.connect(self.RDselect)
        self.RB3.toggled.connect(self.RDselect)
        # here
        self.SuggestTipushButton.clicked.connect(self.TIpushed)
        self.CalT1pushButton.clicked.connect(self.CalT1pushed)
        self.CalStdPMRpushButton.clicked.connect(self.CalPMRpushed)
        PMRcalculater.setCentralWidget(self.widget)

        self.retranslateUi()
        QtCore.QMetaObject.connectSlotsByName(PMRcalculater)

    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.groupBox.setTitle(_translate("PMRcalculater", "Scanner Type"))
        self.RB1_5.setText(_translate("PMRcalculater", "1.5T"))
        self.RB3.setText(_translate("PMRcalculater", "3T"))
        self.groupBox_2.setTitle(_translate("PMRcalculater", "Standard Properties"))
        self.SPInversionTimeL.setText(_translate("PMRcalculater", "Inversion Time"))
        self.SPSegL.setText(_translate("PMRcalculater", "Number of Segments"))
        self.SPEchosL.setText(_translate("PMRcalculater", "Echo Spacing"))
        self.SPFaL.setText(_translate("PMRcalculater", "Flip Angle"))
        self.SPHrL.setText(_translate("PMRcalculater", "Heart Rate"))
        self.SPT1myoL.setText(_translate("PMRcalculater", "T1 myocardium"))
        self.groupBox_3.setTitle(_translate("PMRcalculater", "Step 1: Enter Scanner Properties"))
        self.HrL.setText(_translate("PMRcalculater", "Heart Rate"))
        self.SegL.setText(_translate("PMRcalculater", "Number of Segments"))
        self.FaL.setText(_translate("PMRcalculater", "Flip Angle"))
        self.EchosL.setText(_translate("PMRcalculater", "Echo Spacing"))
        self.TiL.setText(_translate("PMRcalculater", "Inversion Time"))
        self.SuggestTipushButton.setText(_translate("PMRcalculater", "Suggest TI"))
        self.groupBox_4.setTitle(_translate("PMRcalculater", "Step 2: Enter Scanner Measured PMR"))
        self.ScannerPMRL.setText(_translate("PMRcalculater", "Scanner PMR"))
        self.groupBox_5.setTitle(_translate("PMRcalculater", "Step 3: Calculate Plaque T1"))
        self.PT1L.setText(_translate("PMRcalculater", "Plaque T1"))
        self.CalT1pushButton.setText(_translate("PMRcalculater", "Calculate Plaque T1"))
        self.groupBox_6.setTitle(_translate("PMRcalculater", "Step 4: Calculate Standardized PMR"))
        self.CalStdPMRpushButton.setText(_translate("PMRcalculater", "Calculate Standardized PMR"))
        self.StdPMRL.setText(_translate("PMRcalculater", "Standard PMR"))


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    PMRcalculater = QtWidgets.QMainWindow()
    ui = Ui_PMRcalculater()
    ui.setupUi(PMRcalculater)
    PMRcalculater.show()
    sys.exit(app.exec_())
