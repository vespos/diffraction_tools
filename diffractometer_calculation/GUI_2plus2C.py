import sys
import numpy as np
from PyQt5 import QtCore
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QDoubleSpinBox, QLineEdit, QLabel,
    QGroupBox, QHBoxLayout, QVBoxLayout, QGridLayout, QComboBox)
import IPython

import diffAngles as diff


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        """------------------------- CREATE WIDGETS -------------------------"""
        self.initUI()
        self.create_hklGroupBox()
        self.create_NGroupBox()
        self.create_latticeGroupBox()
        self.create_anglesGroupbox()
        self.create_QQGroupbox()
        self.create_OtherGroupbox()


        """----------------------------- LAYOUT -----------------------------"""
        mainLayout = QGridLayout()
        self.setLayout(mainLayout)
        mainLayout.addLayout(self.quitBox, 4, 1)
        mainLayout.addWidget(self.hklGroupBox, 0, 0)
        mainLayout.addWidget(self.NGroupBox, 1, 0)
        mainLayout.addWidget(self.latticeGroupBox, 2, 0)
        mainLayout.addWidget(self.OtherGroupBox, 3, 0)

        mainLayout.addWidget(self.anglesGroupBox, 2, 1)
        mainLayout.addWidget(self.QQGroupBox, 3, 1)


        """------------------------- UPDATE FIELDS -------------------------"""
        self.h.valueChanged.connect(self.update_angles)
        self.k.valueChanged.connect(self.update_angles)
        self.l.valueChanged.connect(self.update_angles)

        self.a.valueChanged.connect(self.update_angles)
        self.b.valueChanged.connect(self.update_angles)
        self.c.valueChanged.connect(self.update_angles)
        self.aa.valueChanged.connect(self.update_angles)
        self.ba.valueChanged.connect(self.update_angles)
        self.ca.valueChanged.connect(self.update_angles)

        self.E.valueChanged.connect(self.update_angles)
        self.alp.valueChanged.connect(self.update_angles)
        

        self.gamma.valueChanged.connect(self.update_hkl)
        self.delta.valueChanged.connect(self.update_hkl)
        self.omega.valueChanged.connect(self.update_hkl)
        self.alpha.valueChanged.connect(self.update_hkl)

        self.show()



    def initUI(self):
        self.setGeometry(300, 500, 100, 150)
        self.setWindowTitle('2+2C diffractometer calculation')

        qbtn = QPushButton('Quit', self)
        qbtn.clicked.connect(QApplication.instance().quit)
        qbtn.resize(qbtn.sizeHint())

        vbox = QVBoxLayout()
        # vbox.addStretch(1)
        vbox.addWidget(qbtn)
        self.quitBox = QHBoxLayout()
        # self.quitBox.addStretch(1)
        self.quitBox.addLayout(vbox)


    def create_hklGroupBox(self):
        self.hklGroupBox = QGroupBox("h,k,l")
        self.h = QDoubleSpinBox(self)
        self.h.setRange(-20,20)
        self.h.setSingleStep(0.1)
        self.k = QDoubleSpinBox(self)
        self.k.setRange(-20,20)
        self.k.setSingleStep(0.1)
        self.l = QDoubleSpinBox(self)
        self.l.setRange(-20,20)
        self.l.setSingleStep(0.1)

        layout = QHBoxLayout()
        layout.addWidget(self.h)
        layout.addWidget(self.k)
        layout.addWidget(self.l)
        self.hklGroupBox.setLayout(layout)


    def create_NGroupBox(self):
        self.NGroupBox = QGroupBox("Surface normal")
        self.Nh = QDoubleSpinBox(self)
        self.Nh.setRange(-20,20)
        self.Nh.setSingleStep(0.1)
        self.Nk = QDoubleSpinBox(self)
        self.Nk.setRange(-20,20)
        self.Nk.setSingleStep(0.1)
        self.Nl = QDoubleSpinBox(self)
        self.Nl.setRange(-20,20)
        self.Nl.setSingleStep(0.1)

        layout = QHBoxLayout()
        layout.addWidget(self.Nh)
        layout.addWidget(self.Nk)
        layout.addWidget(self.Nl)
        self.NGroupBox.setLayout(layout)


    def create_latticeGroupBox(self):
        abcGroupBox = QGroupBox("a,b,c")
        self.a = QDoubleSpinBox(self)
        self.a.setRange(0,100)
        self.a.setSingleStep(0.01)
        self.b = QDoubleSpinBox(self)
        self.b.setRange(0,100)
        self.b.setSingleStep(0.01)
        self.c = QDoubleSpinBox(self)
        self.c.setRange(0,100)
        self.c.setSingleStep(0.01)
        layout = QHBoxLayout()
        layout.addWidget(self.a)
        layout.addWidget(self.b)
        layout.addWidget(self.c)
        abcGroupBox.setLayout(layout)

        latticeAnglesGroupBox = QGroupBox("alpha,beta,gamma")
        self.aa = QDoubleSpinBox(self)
        self.aa.setRange(0,180)
        self.aa.setSingleStep(1)
        self.aa.setValue(90)
        self.ba = QDoubleSpinBox(self)
        self.ba.setRange(0,180)
        self.ba.setSingleStep(1)
        self.ba.setValue(90)
        self.ca = QDoubleSpinBox(self)
        self.ca.setRange(0,180)
        self.ca.setSingleStep(1)
        self.ca.setValue(90)
        layout = QHBoxLayout()
        layout.addWidget(self.aa)
        layout.addWidget(self.ba)
        layout.addWidget(self.ca)
        latticeAnglesGroupBox.setLayout(layout)

        self.latticeGroupBox = QGroupBox("Lattice parameters")
        layout = QVBoxLayout()
        layout.addWidget(abcGroupBox)
        layout.addWidget(latticeAnglesGroupBox)
        self.latticeGroupBox.setLayout(layout)


    def create_OtherGroupbox(self):
        self.OtherGroupBox = QGroupBox("Other parameters")
        layout = QGridLayout()

        E_label = QLabel()
        E_label.setText('Energy (keV)')
        self.E = QDoubleSpinBox(self)
        self.E.setRange(0,25)
        self.E.setSingleStep(0.1)
        self.E.setValue(7)
        layout.addWidget(E_label, 0, 0)
        layout.addWidget(self.E, 0, 1)

        alp_label = QLabel()
        alp_label.setText('Incident angle')
        self.alp = QDoubleSpinBox(self)
        self.alp.setRange(-180,180)
        self.alp.setSingleStep(0.1)
        layout.addWidget(alp_label, 1, 0)
        layout.addWidget(self.alp, 1, 1)

        mode_label = QLabel()
        mode_label.setText('Mode')
        self.modeCombo = QComboBox(self)
        self.modeCombo.addItem("H, fixed incident")
        self.modeCombo.addItem("V, fixed incident")
        self.modeCombo.addItem("H, symmetric")
        self.modeCombo.addItem("V, symmetric")
        layout.addWidget(mode_label, 2, 0)
        layout.addWidget(self.modeCombo, 2, 1)

        self.OtherGroupBox.setLayout(layout)

        
    def create_anglesGroupbox(self):
        self.anglesGroupBox = QGroupBox("Diffraction angles")
        layout = QGridLayout()

        ttheta_label = QLabel()
        ttheta_label.setText('2-theta')
        self.ttheta = QLineEdit("ttheta")
        self.ttheta.setReadOnly(True)
        layout.addWidget(ttheta_label, 0, 0)
        layout.addWidget(self.ttheta, 0, 1)

        gamma_label = QLabel()
        gamma_label.setText('gamma')
        self.gamma = QDoubleSpinBox(self)
        self.gamma.setRange(-360,360)
        self.gamma.setSingleStep(0.5)
        layout.addWidget(gamma_label, 1, 0)
        layout.addWidget(self.gamma, 1, 1)

        delta_label = QLabel()
        delta_label.setText('delta')
        self.delta = QDoubleSpinBox(self)
        self.delta.setRange(-360,360)
        self.delta.setSingleStep(0.5)
        layout.addWidget(delta_label, 2, 0)
        layout.addWidget(self.delta, 2, 1)

        omega_label = QLabel()
        omega_label.setText('omega')
        self.omega = QDoubleSpinBox(self)
        self.omega.setRange(-360,360)
        self.omega.setSingleStep(0.5)
        layout.addWidget(omega_label, 3, 0)
        layout.addWidget(self.omega, 3, 1)

        alpha_label = QLabel()
        alpha_label.setText('alpha')
        self.alpha = QDoubleSpinBox(self)
        self.alpha.setRange(-360,360)
        self.alpha.setSingleStep(0.5)
        layout.addWidget(alpha_label, 4, 0)
        layout.addWidget(self.alpha, 4, 1)

        self.anglesGroupBox.setLayout(layout)


    def create_QQGroupbox(self):
        self.QQGroupBox = QGroupBox("Others")
        layout = QGridLayout()

        QQ_label = QLabel()
        QQ_label.setText('Q')
        self.QQ = QLineEdit("Q")
        self.QQ.setReadOnly(True)
        layout.addWidget(QQ_label, 0, 0)
        layout.addWidget(self.QQ, 0, 1)

        angout_label = QLabel()
        angout_label.setText('angle out')
        self.angout = QLineEdit("angle out")
        self.angout.setReadOnly(True)
        layout.addWidget(angout_label, 1, 0)
        layout.addWidget(self.angout, 1, 1)

        self.QQGroupBox.setLayout(layout)


    def modeComboActivated(self):
        modeTxt = self.modeCombo.currentText()

        if modeTxt == 'H, fixed incident':
            mode = 1
        elif modeTxt == 'V, fixed incident':
            mode = 2

        return mode


    @QtCore.pyqtSlot()
    def update_angles(self):
        hkl = np.array([self.h.value(), self.k.value(), self.l.value()])
        N = np.array([self.Nh.value(), self.Nk.value(), self.Nl.value()])
        a = np.array([self.a.value(), self.b.value(), self.c.value()])
        aa = np.array([self.aa.value(), self.ba.value(), self.ca.value()])

        alp = self.alp.value()
        E = self.E.value()

        mode = self.modeComboActivated()

        delta, gamma, omega, tt, angout, Q = twoPlusTwoC(hkl, a, aa, N, E, mode, alp)

        # IPython.embed()

        for spinbox, value in zip ((self.gamma, self.delta, self.omega, self.alpha), (gamma, delta, omega, alp)):
            spinbox.blockSignals(True)
            spinbox.setValue(value)
            spinbox.blockSignals(False)
        
        for textbox, value in zip ((self.ttheta, self.angout, self.QQ), (tt, angout, Q)):
            string = "{:.2f}".format(value)
            textbox.setText(string)


    @QtCore.pyqtSlot()
    def update_hkl(self):
        delta = self.delta.value()
        gamma = self.gamma.value()
        omega = self.omega.value()
        alpha = self.alpha.value()
        N = np.array([self.Nh.value(), self.Nk.value(), self.Nl.value()])
        a = np.array([self.a.value(), self.b.value(), self.c.value()])
        aa = np.array([self.aa.value(), self.ba.value(), self.ca.value()])

        E = self.E.value()

        mode = self.modeComboActivated()

        hkl, Q, tt = twoPlusTwoC_hkl(E, delta, gamma, omega, alpha, a, aa, N, mode)

        for spinbox, value in zip ((self.h, self.k, self.l), (hkl[0], hkl[1], hkl[2])):
            spinbox.blockSignals(True)
            spinbox.setValue(value)
            spinbox.blockSignals(False)

        for textbox, value in zip ((self.ttheta, self.QQ), (tt, Q)):
            string = "{:.2f}".format(value)
            textbox.setText(string)



def fourC(hkl, a, aa, N, E, *args):
    """ args[0]: incident angle """
    return diff.main_22C(hkl, a, aa, N, E, mode, *args)

def twoPlusTwoC(hkl, a, aa, N, E, mode, *args):
    """ args[0]: incident angle """
    return diff.main_22C(hkl, a, aa, N, E, mode, *args)

def twoPlusTwoC_hkl(E, delta, gamma, omega, alpha, a, aa, N, mode, *args):
    return diff.main_22C_hkl(E, delta, gamma, omega, alpha, a, aa, N, mode, *args)



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())