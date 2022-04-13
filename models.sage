## The Fricke parameterization for X_0(n)

j2_1(t) = (t+256)^3/t^2
j2_2(t) = (t+16)^3/t
j3_1(t) = (t+27)*(t+243)^3/t^3
j3_2(t) = (t+27)*(t+3)^3/t
j4_1(t) = (t^2+256*t+4096)^3/(t^4*(t+16))
j4_2(t) = (t^2+16*t+16)^3/(t*(t+16))
j5_1(t) = (t^2+250*t+3125)^3/t^5
j5_2(t) = (t^2+10*t+5)^3/t
j6_1(t) = (t+12)^3*(t^3+252*t^2+3888*t+15552)^3/(t^6*(t+8)^2*(t+9)^3)
j6_2(t) = (t+6)^3*(t^3+18*t^2+84*t+24)^3/(t*(t+8)^3*(t+9)^2)
j7_1(t) = (t^2+13*t+49)*(t^2+245*t+2401)^3/t^7
j7_2(t) = (t^2+13*t+49)*(t^2+5*t+1)^3/t
j8_1(t) = (t^4 + 240*t^3 +2144*t^2 +3840*t+256)^3/((t-4)^8*t*(t+4)^2)
j8_2(t) = (t^4 - 16*t^2+16)^3/((t^2-16)*t^2)
j9_1(t) = (t+6)^3*(t^3+234*t^2+756*t+2160)^3/((t-3)^8*(t^3-27))
j9_2(t) = t^3*(t^3-24)^3/(t^3-27)
j10_1(t) = (t^6+236*t^5+1440*t^4+1920*t^3+3840*t^2+256*t+256)^3/(t^2*(t-4)^10*(t+1)^5)
j10_2(t) = (t^6-4*t^5+16*t+16)^3/(t^5*(t-4)*(t+1)^2)
j12_1(t) = (t^2+6*t-3)^3*(t^6+234*t^5+747*t^4+540*t^3-729*t^2-486*t-243)^3/((t-3)^12*(t-1)*t^3*(t+1)^4*(t+3)^3)
j12_2(t) = (t^2-3)^3*(t^6-9*t^4+3*t^2-3)^3/(t^4*(t^2-9)*(t^2-1)^3)
j13_1(t) = (t^2+5*t+13)*(t^4+247*t^3+3380*t^2+15379*t+28561)^3/t^13
j13_2(t) = (t^2+5*t+13)*(t^4+7*t^3+20*t^2+19*t+1)^3/t
j16_1(t) = (t^8+240*t^7+2160*t^6+6720*t^5+17504*t^4+26880*t^3+34560*t^2+15360*t+256)^3/((t-2)^16*t*(t+2)^4*(t^2+4))
j16_2(t) = (t^8-16*t^4+16)^3/(t^4*(t^4-16))
j18_1(t) = (t^3+6*t^2+4)^3*(t^9+234*t^8+756*t^7+2172*t^6+1872*t^5+3024*t^4+48*t^3+3744*t^2+64)^3/((t-2)^18*t^2*(t+1)^9*(t^2-t+1)*(t^2+2*t+4)^2)
j18_2(t) = (t^3-2)^3*(t^9-6*t^6-12*t^3-8)^3/(t^9*(t^3-8)*(t^3+1)^2)
j25_1(t) = (t - 1)^-25 * (t^4 + t^3 + 6*t^2 + 6*t + 11)^-1 * (t^10 + 240*t^9 + 2170*t^8 + 8880*t^7 + 34835*t^6 + 83748*t^5 + 206210*t^4 + 313380*t^3 + 503545*t^2 + 424740*t + 375376)^3
j25_2(t) = (t^10+10*t^8+35*t^6-12*t^5+50*t^4-60*t^3+25*t^2-60*t+16)^3/(t^5+5*t^3+5*t-11)

## The Models \mathcal{C}_{n,i}(t,d)

def calC2_1(t,d):
	A=(-27) * (t + 64) * (t + 256) * d^2
	B=(-54) * (t - 512) * (t + 64)^2 * d^3
	return EllipticCurve([A,B])
def calC2_2(t,d):
	A=(-432) * (t + 16) * (t + 64) * d^2
	B=(-3456) * (t - 8) * (t + 64)^2 * d^3
	return EllipticCurve([A,B])
def calC3_1(t,d):
	A=(-3) * (t + 243) * d^2 * (t + 27)^3
	B=(-2) * d^3 * (t + 27)^4 * (t^2 - 486*t - 19683)
	return EllipticCurve([A,B])
def calC3_2(t,d):
	A=(-243) * (t + 3) * d^2 * (t + 27)^3
	B=(-1458) * d^3 * (t + 27)^4 * (t^2 + 18*t - 27)
	return EllipticCurve([A,B])
def calC4_1(t,d):
	A=(-432) * d^2 * (t^2 + 16*t + 256)
	B=(3456) * (t - 16) * (t + 8) * (t + 32) * d^3
	return EllipticCurve([A,B])
def calC4_2(t,d):
	A=(-6912) * d^2 * (t^2 + 16*t + 16)
	B=(221184) * (t + 8) * d^3 * (t^2 + 16*t - 8)
	return EllipticCurve([A,B])
def calC4_3(t,d):
	A=(-432) * d^2 * (t^2 - 224*t + 256)
	B=(3456) * (t - 16) * d^3 * (t^2 + 544*t + 256)
	return EllipticCurve([A,B])
def calC4_4(t,d):
	A=(-27) * d^2 * (t^2 + 256*t + 4096)
	B=(54) * (t + 32) * d^3 * (t^2 - 512*t - 8192)
	return EllipticCurve([A,B])
def calC5_1(t,d):
	A=(-27) * d^2 * (t^2 + 22*t + 125) * (t^2 + 250*t + 3125)
	B=(-54) * d^3 * (t^2 - 500*t - 15625) * (t^2 + 22*t + 125)^2
	return EllipticCurve([A,B])
def calC5_2(t,d):
	A=(-16875) * d^2 * (t^2 + 10*t + 5) * (t^2 + 22*t + 125)
	B=(-843750) * d^3 * (t^2 + 4*t - 1) * (t^2 + 22*t + 125)^2
	return EllipticCurve([A,B])
def calC6_1(t,d):
	A=(-3) * (t + 12) * d^2 * (t^3 + 252*t^2 + 3888*t + 15552)
	B=(-2) * d^3 * (t^2 + 36*t + 216) * (t^4 - 504*t^3 - 13824*t^2 - 124416*t - 373248)
	return EllipticCurve([A,B])
def calC6_2(t,d):
	A=(-48) * (t + 6) * d^2 * (t^3 + 18*t^2 + 324*t + 1944)
	B=(-128) * d^3 * (t^2 + 36*t + 216) * (t^4 - 216*t^2 - 1944*t - 5832)
	return EllipticCurve([A,B])
def calC6_3(t,d):
	A=(-243) * (t + 12) * d^2 * (t^3 + 12*t^2 + 48*t + 192)
	B=(-1458) * d^3 * (t^2 + 12*t + 24) * (t^4 + 24*t^3 + 192*t^2 - 4608)
	return EllipticCurve([A,B])
def calC6_4(t,d):
	A=(-3888) * (t + 6) * d^2 * (t^3 + 18*t^2 + 84*t + 24)
	B=(-93312) * d^3 * (t^2 + 12*t + 24) * (t^4 + 24*t^3 + 192*t^2 + 504*t - 72)
	return EllipticCurve([A,B])
def calC7_1(t,d):
	A=(-27) * d^2 * (t^2 + 13*t + 49) * (t^2 + 245*t + 2401)
	B=(-54) * d^3 * (t^2 + 13*t + 49) * (t^4 - 490*t^3 - 21609*t^2 - 235298*t - 823543)
	return EllipticCurve([A,B])
def calC7_2(t,d):
	A=(-64827) * d^2 * (t^2 + 5*t + 1) * (t^2 + 13*t + 49)
	B=(-6353046) * d^3 * (t^2 + 13*t + 49) * (t^4 + 14*t^3 + 63*t^2 + 70*t - 7)
	return EllipticCurve([A,B])
def calC8_1(t,d):
	A=(-432) * d^2 * (t^4 + 224*t^2 + 256)
	B=(-3456) * d^3 * (t^2 + 16) * (t^2 - 24*t + 16) * (t^2 + 24*t + 16)
	return EllipticCurve([A,B])
def calC8_2(t,d):
	A=(-432) * d^2 * (t^4 - 240*t^3 + 2144*t^2 - 3840*t + 256)
	B=(-3456) * d^3 * (t^2 - 24*t + 16) * (t^4 + 528*t^3 - 4000*t^2 + 8448*t + 256)
	return EllipticCurve([A,B])
def calC8_3(t,d):
	A=(-27) * d^2 * (t^4 + 240*t^3 + 2144*t^2 + 3840*t + 256)
	B=(-54) * d^3 * (t^2 + 24*t + 16) * (t^4 - 528*t^3 - 4000*t^2 - 8448*t + 256)
	return EllipticCurve([A,B])
def calC8_4(t,d):
	A=(-6912) * d^2 * (t^4 - 16*t^2 + 256)
	B=(-221184) * d^3 * (t^2 - 32) * (t^2 - 8) * (t^2 + 16)
	return EllipticCurve([A,B])
def calC8_5(t,d):
	A=(-6912) * d^2 * (t^4 - 256*t^2 + 4096)
	B=(-221184) * d^3 * (t^2 - 32) * (t^4 + 512*t^2 - 8192)
	return EllipticCurve([A,B])
def calC8_6(t,d):
	A=(-110592) * d^2 * (t^4 - 16*t^2 + 16)
	B=(-14155776) * d^3 * (t^2 - 8) * (t^4 - 16*t^2 - 8)
	return EllipticCurve([A,B])
def calC9_1(t,d):
	A=(-3) * (t + 6) * d^2 * (t^3 + 234*t^2 + 756*t + 2160)
	B=(-2) * d^3 * (t^6 - 504*t^5 - 16632*t^4 - 123012*t^3 - 517104*t^2 - 1143072*t - 1475496)
	return EllipticCurve([A,B])
def calC9_2(t,d):
	A=(-243) * t * (t + 6) * d^2 * (t^2 - 6*t + 36)
	B=(-1458) * d^3 * (t^2 - 6*t - 18) * (t^4 + 6*t^3 + 54*t^2 - 108*t + 324)
	return EllipticCurve([A,B])
def calC9_3(t,d):
	A=(-19683) * t * d^2 * (t^3 - 24)
	B=(-1062882) * d^3 * (t^6 - 36*t^3 + 216)
	return EllipticCurve([A,B])
def calC10_1(t,d):
	A=(-27) * d^2 * (t^2 + 4) * (t^6 + 236*t^5 + 1440*t^4 + 1920*t^3 + 3840*t^2 + 256*t + 256)
	B=(-54) * d^3 * (t^2 + 4*t + 8) * (t^2 + 22*t - 4) * (t^2 + 4)^2 * (t^4 - 536*t^3 - 264*t^2 - 416*t - 64)
	return EllipticCurve([A,B])
def calC10_2(t,d):
	A=(-432) * d^2 * (t^2 + 4) * (t^6 - 4*t^5 + 240*t^4 - 480*t^3 + 1440*t^2 - 944*t + 16)
	B=(-3456) * d^3 * (t^2 - 2*t + 2) * (t^2 + 22*t - 4) * (t^2 + 4)^2 * (t^4 - 26*t^3 + 66*t^2 - 536*t - 4)
	return EllipticCurve([A,B])
def calC10_3(t,d):
	A=(-16875) * d^2 * (t^2 + 4) * (t^6 - 4*t^5 + 256*t + 256)
	B=(-843750) * d^3 * (t^2 - 2*t - 4) * (t^2 + 4*t + 8) * (t^2 + 4)^2 * (t^4 - 8*t^3 + 24*t^2 - 32*t - 64)
	return EllipticCurve([A,B])
def calC10_4(t,d):
	A=(-270000) * d^2 * (t^2 + 4) * (t^6 - 4*t^5 + 16*t + 16)
	B=(-54000000) * d^3 * (t^2 - 2*t - 4) * (t^2 - 2*t + 2) * (t^2 + 4)^2 * (t^4 - 2*t^3 - 6*t^2 - 8*t - 4)
	return EllipticCurve([A,B])
def calC11_1(t,d):
	A=(-1149984) * d^2
	B=(-487018224) * d^3
	return EllipticCurve([A,B])
def calC11_2(t,d):
	A=(-9504) * d^2
	B=(365904) * d^3
	return EllipticCurve([A,B])
def calC11_3(t,d):
	A=(-395307) * d^2
	B=(373960422) * d^3
	return EllipticCurve([A,B])
def calC11_4(t,d):
	A=(-38907) * d^2
	B=(-2953962)  * d^3
	return EllipticCurve([A,B])
def calC12_1(t,d):
	A=(-48) * d^2 * (t^2 + 3) * (t^6 + 225*t^4 - 405*t^2 + 243)
	B=(-128) * d^3 * (t^4 + 18*t^2 - 27) * (t^4 - 24*t^3 + 18*t^2 - 27) * (t^4 + 24*t^3 + 18*t^2 - 27)
	return EllipticCurve([A,B])
def calC12_2(t,d):
	A=(-3888) * d^2 * (t^2 + 3) * (t^6 - 15*t^4 + 75*t^2 + 3)
	B=(-93312) * d^3 * (t^4 - 6*t^2 - 3) * (t^4 - 6*t^2 - 24*t - 3) * (t^4 - 6*t^2 + 24*t - 3)
	return EllipticCurve([A,B])
def calC12_3(t,d):
	A=(-768) * d^2 * (t^2 - 3) * (t^6 - 9*t^4 + 243*t^2 - 243)
	B=(-8192) * d^3 * (t^4 + 18*t^2 - 27) * (t^8 - 36*t^6 + 270*t^4 - 972*t^2 + 729)
	return EllipticCurve([A,B])
def calC12_4(t,d):
	A=(-62208) * d^2 * (t^2 - 3) * (t^6 - 9*t^4 + 3*t^2 - 3)
	B=(-5971968) * d^3 * (t^4 - 6*t^2 - 3) * (t^8 - 12*t^6 + 30*t^4 - 36*t^2 + 9)
	return EllipticCurve([A,B])
def calC12_5(t,d):
	A=(-3) * d^2 * (t^2 + 6*t - 3) * (t^6 + 234*t^5 + 747*t^4 + 540*t^3 - 729*t^2 - 486*t - 243)
	B=(-2) * d^3 * (t^4 + 24*t^3 + 18*t^2 - 27) * (t^8 - 528*t^7 - 3996*t^6 - 9504*t^5 + 270*t^4 + 14256*t^3 - 972*t^2 + 729)
	return EllipticCurve([A,B])
def calC12_6(t,d):
	A=(-243) * d^2 * (t^2 + 6*t - 3) * (t^6 - 6*t^5 + 27*t^4 + 60*t^3 - 249*t^2 + 234*t - 3)
	B=(-1458) * d^3 * (t^4 - 6*t^2 + 24*t - 3) * (t^8 - 12*t^6 - 528*t^5 + 30*t^4 + 3168*t^3 - 3996*t^2 + 1584*t + 9)
	return EllipticCurve([A,B])
def calC12_7(t,d):
	A=(-48) * d^2 * (t^2 - 6*t - 3) * (t^6 - 234*t^5 + 747*t^4 - 540*t^3 - 729*t^2 + 486*t - 243)
	B=(-128) * d^3 * (t^4 - 24*t^3 + 18*t^2 - 27) * (t^8 + 528*t^7 - 3996*t^6 + 9504*t^5 + 270*t^4 - 14256*t^3 - 972*t^2 + 729)
	return EllipticCurve([A,B])
def calC12_8(t,d):
	A=(-3888) * d^2 * (t^2 - 6*t - 3) * (t^6 + 6*t^5 + 27*t^4 - 60*t^3 - 249*t^2 - 234*t - 3)
	B=(-93312) * d^3 * (t^4 - 6*t^2 - 24*t - 3) * (t^8 - 12*t^6 + 528*t^5 + 30*t^4 - 3168*t^3 - 3996*t^2 - 1584*t + 9)
	return EllipticCurve([A,B])
def calC13_1(t,d):
	A=(-27) * d^2 * (t^2 + 5*t + 13) * (t^2 + 6*t + 13) * (t^4 + 247*t^3 + 3380*t^2 + 15379*t + 28561)
	B=(-54) * d^3 * (t^2 + 5*t + 13) * (t^2 + 6*t + 13)^2 * (t^6 - 494*t^5 - 20618*t^4 - 237276*t^3 - 1313806*t^2 - 3712930*t - 4826809)
	return EllipticCurve([A,B])
def calC13_2(t,d):
	A=(-771147) * d^2 * (t^2 + 5*t + 13) * (t^2 + 6*t + 13) * (t^4 + 7*t^3 + 20*t^2 + 19*t + 1)
	B=(-260647686) * d^3 * (t^2 + 5*t + 13) * (t^2 + 6*t + 13)^2 * (t^6 + 10*t^5 + 46*t^4 + 108*t^3 + 122*t^2 + 38*t - 1)
	return EllipticCurve([A,B])
def calC14_1(d):
	A=(-2361555) * d^2
	B=(1396762542) * d^3
	return EllipticCurve([A,B])
def calC14_2(d):
	A=(-138915) * d^2
	B=(24504606) * d^3
	return EllipticCurve([A,B])
def calC14_3(d):
	A=(-48195) * d^2
	B=(-4072194) * d^3
	return EllipticCurve([A,B])
def calC14_4(d):
	A=(-2835) * d^2
	B=(-71442) * d^3
	return EllipticCurve([A,B])
def calC15_1(d):
	A=(-162675) * d^2
	B=(-25254450) * d^3
	return EllipticCurve([A,B])
def calC15_2(d):
	A=(-675) * d^2
	B=(-79650) * d^3
	return EllipticCurve([A,B])
def calC15_3(d):
	A=(712125) * d^2
	B=(-104861250) * d^3
	return EllipticCurve([A,B])
def calC15_4(d):
	A=(-97875) * d^2
	B=(14208750) * d^3
	return EllipticCurve([A,B])
def calC16_1(t,d):
	A=(-432) * d^2 * (t^8 + 240*t^6 + 2144*t^4 + 3840*t^2 + 256)
	B=(-3456) * d^3 * (t^4 + 24*t^2 + 16) * (t^4 - 24*t^3 + 24*t^2 - 96*t + 16) * (t^4 + 24*t^3 + 24*t^2 + 96*t + 16)
	return EllipticCurve([A,B])
def calC16_2(t,d):
	A=(-27) * d^2 * (t^8 + 240*t^7 + 2160*t^6 + 6720*t^5 + 17504*t^4 + 26880*t^3 + 34560*t^2 + 15360*t + 256)
	B=(-54) * d^3 * (t^4 + 24*t^3 + 24*t^2 + 96*t + 16) * (t^8 - 528*t^7 - 3984*t^6 - 14784*t^5 - 31648*t^4 - 59136*t^3 - 63744*t^2 - 33792*t + 256)
	return EllipticCurve([A,B])
def calC16_3(t,d):
	A=(-432) * d^2 * (t^8 - 240*t^7 + 2160*t^6 - 6720*t^5 + 17504*t^4 - 26880*t^3 + 34560*t^2 - 15360*t + 256)
	B=(-3456) * d^3 * (t^4 - 24*t^3 + 24*t^2 - 96*t + 16) * (t^8 + 528*t^7 - 3984*t^6 + 14784*t^5 - 31648*t^4 + 59136*t^3 - 63744*t^2 + 33792*t + 256)
	return EllipticCurve([A,B])
def calC16_4(t,d):
	A=(-6912) * d^2 * (t^4 - 4*t^3 + 8*t^2 + 16*t + 16) * (t^4 + 4*t^3 + 8*t^2 - 16*t + 16)
	B=(-221184) * d^3 * (t^2 - 4*t - 4) * (t^2 + 4*t - 4) * (t^4 + 16) * (t^4 + 24*t^2 + 16)
	return EllipticCurve([A,B])
def calC16_5(t,d):
	A=(-6912) * d^2 * (t^4 - 16*t^3 + 8*t^2 + 64*t + 16) * (t^4 + 16*t^3 + 8*t^2 - 64*t + 16)
	B=(-221184) * d^3 * (t^2 - 4*t - 4) * (t^2 + 4*t - 4) * (t^8 + 528*t^6 - 4000*t^4 + 8448*t^2 + 256)
	return EllipticCurve([A,B])
def calC16_6(t,d):
	A=(-110592) * d^2 * (t^8 - 16*t^4 + 256)
	B=(-14155776) * d^3 * (t^4 - 32) * (t^4 - 8) * (t^4 + 16)
	return EllipticCurve([A,B])
def calC16_7(t,d):
	A=(-110592) * d^2 * (t^8 - 256*t^4 + 4096)
	B=(-14155776) * d^3 * (t^4 - 32) * (t^8 + 512*t^4 - 8192)
	return EllipticCurve([A,B])
def calC16_8(t,d):
	A=(-1769472) * d^2 * (t^8 - 16*t^4 + 16)
	B=(-905969664) * d^3 * (t^4 - 8) * (t^8 - 16*t^4 - 8)
	return EllipticCurve([A,B])
def calC17_1(d):
	A=(-247394115) * d^2
	B=(-1679010134850) * d^3
	return EllipticCurve([A,B])
def calC17_2(d):
	A=(-3940515) * d^2
	B=(3010787550) * d^3
	return EllipticCurve([A,B])
def calC18_1(t,d):
	A=(-3) * d^2 * (t^3 + 6*t^2 + 4) * (t^9 + 234*t^8 + 756*t^7 + 2172*t^6 + 1872*t^5 + 3024*t^4 + 48*t^3 + 3744*t^2 + 64)
	B=(-2) * d^3 * (t^6 + 24*t^5 + 24*t^4 + 92*t^3 - 48*t^2 + 96*t - 8) * (t^12 - 528*t^11 - 3984*t^10 - 14792*t^9 - 27936*t^8 - 42624*t^7 - 37632*t^6 - 52992*t^5 - 25344*t^4 - 43520*t^3 - 6144*t^2 - 6144*t - 512)
	return EllipticCurve([A,B])
def calC18_2(t,d):
	A=(-48) * d^2 * (t^3 + 6*t - 2) * (t^9 + 234*t^7 - 6*t^6 + 756*t^5 - 936*t^4 + 2172*t^3 - 1512*t^2 + 936*t - 8)
	B=(-128) * d^3 * (t^6 + 24*t^5 + 24*t^4 + 92*t^3 - 48*t^2 + 96*t - 8) * (t^12 - 24*t^11 + 48*t^10 - 680*t^9 + 792*t^8 - 3312*t^7 + 4704*t^6 - 10656*t^5 + 13968*t^4 - 14792*t^3 + 7968*t^2 - 2112*t - 8)
	return EllipticCurve([A,B])
def calC18_3(t,d):
	A=(-243) * d^2 * (t^3 + 4) * (t^3 + 6*t^2 + 4) * (t^6 - 6*t^5 + 36*t^4 + 8*t^3 - 24*t^2 + 16)
	B=(-1458) * d^3 * (t^2 + 2*t - 2) * (t^4 - 8*t^3 - 8*t - 8) * (t^4 - 2*t^3 + 6*t^2 + 4*t + 4) * (t^8 + 8*t^7 + 64*t^6 - 16*t^5 - 56*t^4 + 128*t^3 + 64*t^2 - 64*t + 64)
	return EllipticCurve([A,B])
def calC18_4(t,d):
	A=(-3888) * d^2 * (t^3 - 2) * (t^3 + 6*t - 2) * (t^6 - 6*t^4 - 4*t^3 + 36*t^2 + 12*t + 4)
	B=(-93312) * d^3 * (t^2 + 2*t - 2) * (t^4 - 2*t^3 - 8*t - 2) * (t^4 - 2*t^3 + 6*t^2 + 4*t + 4) * (t^8 + 2*t^7 + 4*t^6 - 16*t^5 - 14*t^4 + 8*t^3 + 64*t^2 - 16*t + 4)
	return EllipticCurve([A,B])
def calC18_5(t,d):
	A=(-19683) * d^2 * (t^3 + 4) * (t^9 - 12*t^6 + 48*t^3 + 64)
	B=(-1062882) * d^3 * (t^6 - 4*t^3 - 8) * (t^12 - 8*t^9 - 512*t^3 - 512)
	return EllipticCurve([A,B])
def calC18_6(t,d):
	A=(-314928) * d^2 * (t^3 - 2) * (t^9 - 6*t^6 - 12*t^3 - 8)
	B=(-68024448) * d^3 * (t^6 - 4*t^3 - 8) * (t^12 - 8*t^9 - 8*t^3 - 8)
	return EllipticCurve([A,B])
def calC19_1(d):
	A=(-219488) * d^2
	B=(-39617584) * d^3
	return EllipticCurve([A,B])
def calC19_2(d):
	A=(-608) * d^2
	B=(5776) * d^3
	return EllipticCurve([A,B])
def calC21_1(d):
	A=(-1396035) * d^2
	B=(634881726) * d^3
	return EllipticCurve([A,B])
def calC21_2(d):
	A=(-1104435) * d^2
	B=(907504398) * d^3
	return EllipticCurve([A,B])
def calC21_3(d):
	A=(3645) * d^2
	B=(-13122) * d^3
	return EllipticCurve([A,B])
def calC21_4(d):
	A=(-54675) * d^2
	B=(-5156946) * d^3
	return EllipticCurve([A,B])
def calC25_1(t,d):
	A=(-27) * d^2 * (t^2 + 4) * (t^10 + 240*t^9 + 2170*t^8 + 8880*t^7 + 34835*t^6 + 83748*t^5 + 206210*t^4 + 313380*t^3 + 503545*t^2 + 424740*t + 375376)
	B=(-54) * d^3 * (t^2 + 4)^2 * (t^4 + 6*t^3 + 21*t^2 + 36*t + 61) * (t^10 - 510*t^9 - 13580*t^8 - 36870*t^7 - 190915*t^6 - 393252*t^5 - 1068040*t^4 - 1508370*t^3 - 2581955*t^2 - 2087010*t - 1885124)
	return EllipticCurve([A,B])
def calC25_2(t,d):
	A=(-16875) * d^2 * (t^2 + 4) * (t^2 + 3*t + 1) * (t^4 - 4*t^3 + 11*t^2 - 14*t + 31) * (t^4 + t^3 + 11*t^2 - 4*t + 16)
	B=(-843750) * d^3 * (t^2 - 2*t - 4) * (t^2 + 4)^2 * (t^4 + 3*t^2 + 1) * (t^4 - 4*t^3 + 21*t^2 - 34*t + 41) * (t^4 + 6*t^3 + 21*t^2 + 36*t + 61)
	return EllipticCurve([A,B])
def calC25_3(t,d):
	A=(-10546875) * d^2 * (t^2 + 4) * (t^10 + 10*t^8 + 35*t^6 - 12*t^5 + 50*t^4 - 60*t^3 + 25*t^2 - 60*t + 16)
	B=(-13183593750) * d^3 * (t^2 + 4)^2 * (t^4 + 3*t^2 + 1) * (t^10 + 10*t^8 + 35*t^6 - 18*t^5 + 50*t^4 - 90*t^3 + 25*t^2 - 90*t + 76)
	return EllipticCurve([A,B])
def calC27_1(d):
	A=(-4320) * d^2
	B=(-109296) * d^3
	return EllipticCurve([A,B])
def calC27_2(d):
	A=0
	B=- 432*d^3
	return EllipticCurve([A,B])
def calC27_3(d):
	A=0
	B=16*d^3
	return EllipticCurve([A,B])
def calC27_4(d):
	A=(-480) * d^2
	B=(4048) * d^3
	return EllipticCurve([A,B])
def calC37_1(d):
	A=(-269675595) * d^2
	B=(-1704553285050) * d^3
	return EllipticCurve([A,B])
def calC37_2(d):
	A=(-10395) * d^2
	B=(444150) * d^3
	return EllipticCurve([A,B])
def calC43_1(d):
	A=(-25442240) * d^2
	B=(-49394836848) * d^3
	return EllipticCurve([A,B])
def calC43_2(d):
	A=(-13760) * d^2
	B=(621264) * d^3
	return EllipticCurve([A,B])
def calC67_1(d):
	A=(-529342880) * d^2
	B=(-4687634371504) * d^3
	return EllipticCurve([A,B])
def calC67_2(d):
	A=(-117920) * d^2
	B=(15585808) * d^3
	return EllipticCurve([A,B])
def calC163_1(d):
	A=(-924354639680) * d^2
	B=(-342062961763303088) * d^3
	return EllipticCurve([A,B])
def calC163_2(d):
	A=(-34790720) * d^2
	B=(78984748304) * d^3
	return EllipticCurve([A,B])

## The j-invariant of \mathcal{C}_{n,i}(t,d), which we define as Jni(t):

J21(t)=factor(calC2_1(t,d).j_invariant())
J22(t)=factor(calC2_2(t,d).j_invariant())
J31(t)=factor(calC3_1(t,d).j_invariant())
J32(t)=factor(calC3_2(t,d).j_invariant())
J41(t)=factor(calC4_1(t,d).j_invariant())
J42(t)=factor(calC4_2(t,d).j_invariant())
J43(t)=factor(calC4_3(t,d).j_invariant())
J44(t)=factor(calC4_4(t,d).j_invariant())
J51(t)=factor(calC5_1(t,d).j_invariant())
J52(t)=factor(calC5_2(t,d).j_invariant())
J61(t)=factor(calC6_1(t,d).j_invariant())
J62(t)=factor(calC6_2(t,d).j_invariant())
J63(t)=factor(calC6_3(t,d).j_invariant())
J64(t)=factor(calC6_4(t,d).j_invariant())
J71(t)=factor(calC7_1(t,d).j_invariant())
J72(t)=factor(calC7_2(t,d).j_invariant())
J81(t)=factor(calC8_1(t,d).j_invariant())
J82(t)=factor(calC8_2(t,d).j_invariant())
J83(t)=factor(calC8_3(t,d).j_invariant())
J84(t)=factor(calC8_4(t,d).j_invariant())
J85(t)=factor(calC8_5(t,d).j_invariant())
J86(t)=factor(calC8_6(t,d).j_invariant())
J91(t)=factor(calC9_1(t,d).j_invariant())
J92(t)=factor(calC9_2(t,d).j_invariant())
J93(t)=factor(calC9_3(t,d).j_invariant())
J101(t)=factor(calC10_1(t,d).j_invariant())
J102(t)=factor(calC10_2(t,d).j_invariant())
J103(t)=factor(calC10_3(t,d).j_invariant())
J104(t)=factor(calC10_4(t,d).j_invariant())
J121(t)=factor(calC12_1(t,d).j_invariant())
J122(t)=factor(calC12_2(t,d).j_invariant())
J123(t)=factor(calC12_3(t,d).j_invariant())
J124(t)=factor(calC12_4(t,d).j_invariant())
J125(t)=factor(calC12_5(t,d).j_invariant())
J126(t)=factor(calC12_6(t,d).j_invariant())
J127(t)=factor(calC12_7(t,d).j_invariant())
J128(t)=factor(calC12_8(t,d).j_invariant())
J131(t)=factor(calC13_1(t,d).j_invariant())
J132(t)=factor(calC13_2(t,d).j_invariant())
J161(t)=factor(calC16_1(t,d).j_invariant())
J162(t)=factor(calC16_2(t,d).j_invariant())
J163(t)=factor(calC16_3(t,d).j_invariant())
J164(t)=factor(calC16_4(t,d).j_invariant())
J165(t)=factor(calC16_5(t,d).j_invariant())
J166(t)=factor(calC16_6(t,d).j_invariant())
J167(t)=factor(calC16_7(t,d).j_invariant())
J168(t)=factor(calC16_8(t,d).j_invariant())
J181(t)=factor(calC18_1(t,d).j_invariant())
J182(t)=factor(calC18_2(t,d).j_invariant())
J183(t)=factor(calC18_3(t,d).j_invariant())
J184(t)=factor(calC18_4(t,d).j_invariant())
J185(t)=factor(calC18_5(t,d).j_invariant())
J186(t)=factor(calC18_6(t,d).j_invariant())
J251(t)=factor(calC25_1(t,d).j_invariant())
J252(t)=factor(calC25_2(t,d).j_invariant())
J253(t)=factor(calC25_3(t,d).j_invariant())


## The models F_n=F_n(a,b) for n=6,8,9:
def F4(a,b):
    return EllipticCurve([0, -16*a + b, 0, -16*a*b, 0])
def F6(a,b):
    return EllipticCurve([36*a + 5*b, 18*a*b + 2*b^2, 648*a^2*b + 153*a*b^2 + 9*b^3, 0, 0])
def F8(a,b):
    return EllipticCurve([4*b, -16*a^2 + b^2, -64*a^2*b + 4*b^3, 0, 0])
def F9(a,b):
    return EllipticCurve([18*a + 3*b, 0, -27*a^3 + 27*a^2*b - 9*a*b^2 + b^3, 0, 0])

## The quantities mun and nun for n=6,8,9
mu4=65536*a^4*r - 4096*a^3*b*r - 256*a^2*b^2*r - 4294967296*a^2*b^2*s - 268435456*a*b^3*s + 16777216*b^4*s
nu4=r + 16777216*s
mu6=11019960576*a^11*r - 3673320192*a^10*b*r + 816293376*a^9*b^2*r - 151637832*a^8*b^3*r + 25450119*a^7*b^4*r - 4002210*a^6*b^5*r - 283029961293*a^5*b^6*r - 130118146632*a^4*b^7*r - 22390284303*a^3*b^8*r - 1709554734*a^2*b^9*r - 48876939*a*b^10*r + 6190069752180534021193728*a^5*b^6*s + 2763423996509166973747200*a^4*b^7*s + 461850029046207998853120*a^3*b^8*s + 34259681671399870562304*a^2*b^9*s + 951871709474295644160*a*b^10*s + 2056589122535424*b^11*s
nu6=13976817921*a^3*r + 2931378687*a^2*b*r + 154461213*a*b^2*r + 603419*b^3*r - 305682456897804149194752*a^3*s - 60044768319211529306112*a^2*b*s - 3019928372956977168384*a*b^2*s - 11758205543242530816*b^3*s
mu8=4294967296*a^10*r - 910163968*a^8*b^2*r + 80412672*a^6*b^4*r - 3252224*a^4*b^6*r + 49920*a^2*b^8*r + 54887620458577920*a^8*b^2*s - 13968195719266304*a^6*b^4*s + 1349100767281152*a^4*b^6*s - 59648505806848*a^2*b^8*s + 1099511627776*b^10*s
nu8=-43456*a^2*r - 195*b^2*r - 837518622720*a^2*s - 729070698496*b^2*s
mu9=33352704656118*a^11*r - 12878513915436*a^10*b*r - 4799439950913*a^9*b^2*r + 8548708907664*a^8*b^3*r - 5614961320287*a^7*b^4*r + 2352593883816*a^6*b^5*r - 677531965770*a^5*b^6*r + 133903606464*a^4*b^7*r - 17399973204*a^3*b^8*r + 1341446508*a^2*b^9*r - 46555317*a*b^10*r - 74851166539560535695360*a^11*s + 199936460410954741831680*a^10*b*s - 233760216060366317379072*a^9*b^2*s + 159056456590525171986432*a^8*b^3*s - 72778322377607129501184*a^7*b^4*s + 26187310755391335977664*a^6*b^5*s - 8732547237171576699648*a^5*b^6*s + 2703704233965504337440*a^4*b^7*s - 660493204341440551416*a^3*b^8*s + 108897827301899438256*a^2*b^9*s - 10511342580558830544*a*b^10*s + 450283905890997363*b^11*s
nu9=-1256740569*a^3*r + 580989240*a^2*b*r + 135174780*a*b^2*r + 574757*b^3*r - 67605924155247820800*a^3*s - 34628644525338132480*a^2*b*s - 11112986675198361600*a*b^2*s - 1204404874484944896*b^3*s

