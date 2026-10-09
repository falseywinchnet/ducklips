#!/usr/bin/env python3
"""Small independent contracts for energy-coordinate components and projection."""
import unittest
import numpy as np
from passive_surface import certify,first_order,cascade,delay_line,FlowModel


class EnergyContracts(unittest.TestCase):
    def test_passive_fractional_delay_and_allpass(self):
        for fraction in [0.,.01,.3,.9,1.]:
            fir=first_order([1-fraction,fraction],[1,0])
            self.assertLessEqual(certify(fir),1+1e-14)
            pole=.27
            root=np.sqrt(1-pole*pole)
            ap=(np.array([[-pole]]),np.array([[root]]),np.array([[root]]),pole)
            device=cascade(delay_line(3),cascade(fir,ap))
            self.assertLessEqual(certify(device),1+1e-14)
            a,b,c,d=device
            for omega in [.01,.5,2.9]:
                z=np.exp(1j*omega)
                expected=z**-3*((1-fraction)+fraction/z)*(pole+1/z)/(1+pole/z)
                actual=(c@np.linalg.solve(z*np.eye(len(a))-a,b)).item()+d
                self.assertLess(abs(actual-expected),1e-13)

    def test_active_filter_rejected(self):
        with self.assertRaises(ValueError):
            certify(first_order([1.1,0],[1,0]))
        with self.assertRaises(ValueError):delay_line(0)

    def test_wave_port_closure_has_correct_units(self):
        # Matched anechoic wave port: b=0, so p=U/Y. This catches the factor of
        # two/reference impedance error that a self-comparison would miss.
        flow=FlowModel(np.zeros((1,1)),np.zeros((1,1)),np.zeros((1,1)),np.zeros((1,1)),
                       np.zeros((1,1)),np.zeros((1,1)),np.array([.5]),48000)
        self.assertEqual(flow.D.item(),4.)
        self.assertEqual(flow.response(100)[0],4.)
        self.assertEqual(flow.responses(100)[0,0],4.)

if __name__=='__main__':unittest.main()
