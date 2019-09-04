'''
Ru Meng, Yanwei Zhou, Jianting Ning, Kaitai Liang, Jinguang Han, Willy Susilo
 
| From: "An Efficient Key-Policy Attribute-Based Searchable Encryption in Prime-Order Groups"
| Published in: 
| Available from: 
| Notes: Security Assumption:

* type:           SE supporting monotonic search query
* setting:        Pairing

:Authors:    Yi-Fan Tseng
:Date:            08/30/2019
'''
from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from charm.toolbox.secretutil import SecretUtil
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.hash_module import Waters, Hash
import time

debug = False
class MZNLHS17(ABEnc):
    
    def __init__(self, groupObj):
        ABEnc.__init__(self)
        global util, group
        util = SecretUtil(groupObj, debug)        
        group = groupObj

    def setup(self):
        g, g_hat = group.random(G1), group.random(G2)
        x1, x2, x3 = group.random(), group.random(), group.random()
        u = g ** x1
        u_hat = g_hat ** x1
        h = g ** x2
        h_hat = g_hat ** x2
        delta = g ** x3
        delta_hat = g_hat ** x3
        #u, h, delta = group.random(G1), group.random(G1), group.random(G1)
        #u_hat, h_hat, delta_hat = group.random(G2), group.random(G2), group.random(G2)
        alpha, d1, d2, d3, d4 = group.random(), group.random(), group.random(), group.random(), group.random()
        g1 = g**d1
        g2 = g**d2
        g3 = g**d3
        g4 = g**d4
        pk = {"g":g, "u":u, "h":h, "delta":delta, "g1":g1, "g2":g2, "g3":g3, "g4":g4, "Omega":pair(g, g_hat)**alpha}
        msk = pk.copy()
        msk["alpha"] = alpha
        msk["g_hat"] = g_hat
        msk["u_hat"] = u_hat
        msk["h_hat"] = h_hat
        msk["delta_hat"] = delta_hat
        msk["d1"] = d1
        msk["d2"] = d2
        msk["d3"] = d3
        msk["d4"] = d4

        return (msk, pk)

    def s_keygen(self, pk):
        kappa = group.random()
        (pk_s, sk_s) = (pk["g"]**kappa, kappa)
        return (pk_s, sk_s)

    def encrypt(self, pk, Keyword):
        H = Hash(group)
        mu, s, s1, s2 = group.random(), group.random(), group.random(), group.random()

        C = pk["Omega"]**mu
        C1 = pk["g"]**mu
        C2 = ( pk["h"]  )
        for name, value in Keyword.items():
            C2 *= ( pk["u"] ** H.hashToZr(name + value) )
        C2 = C2 ** s
        C2 *= pk["delta"] ** -mu
        E1 = pk["g1"] ** (s - s1)
        E2 = pk["g2"] ** s1
        E3 = pk["g3"] ** (s - s2)
        E4 = pk["g4"] ** s2        
        
        return {"C":C, "C1":C1, "C2":C2, "E1":E1, "E2":E2, "E3":E3, "E4":E4, "mu":mu}
    
    def keygen(self, msk, pk_s, Pol, Pol_M):
        H = Hash(group)
        SK = dict()
        s, y2, y3 = msk["alpha"], group.random(), group.random()
        r, r_prime = group.random(), group.random()

        D = msk["g"] ** r
        D_hat = msk["g_hat"] ** r_prime
        SK["D"] = D
        SK["D_hat"] = D_hat

        X = pair(pk_s, D_hat) ** r
        X = H.hashToZr(X)
        
        for name, row in Pol_M.items():
            lamb = row[0]*s + row[1]*y2 + row[2]*y3
            t1, t2 = group.random(), group.random()
            d1 = msk["d1"]
            d2 = msk["d2"]
            d3 = msk["d3"]
            d4 = msk["d4"]

            t = d1*d2*t1 + d3*d4*t2
            h = H.hashToZr(name + Pol[name])
            Y = (msk["u_hat"]**h) * msk["h_hat"]

            D = ( msk["g_hat"] ** lamb ) * ( msk["delta_hat"] **  t)
            R = msk["g_hat"]** (t+X)
            #R = msk["g_hat"]** t
            T1 = Y ** (-d2*t1)
            T2 = Y ** (-d1*t1)
            T3 = Y ** (-d4*t2)
            T4 = Y ** (-d3*t2)
            
            Q = dict()
            tmp = Pol.copy()
            del tmp[name]
            for Q_name, Q_value in tmp.items():
                hQ = H.hashToZr(Q_name + Q_value)
                Q_tmp = msk["u_hat"] ** (-hQ)
                Q[Q_name] = {"Q1":Q_tmp ** (d2*t1), "Q2":Q_tmp ** (d1*t1), "Q3":Q_tmp ** (d4*t2), "Q4":Q_tmp ** (d3*t2)}

            SK[name] = {"D":D, "R":R, "T1":T1, "T2":T2, "T3":T3, "T4":T4 ,"Q":Q}
            SK["g_hat"] = msk["g_hat"]

        return SK
    def decrypt(self, pk, sk_s, SK, CT, Delta):
        H = Hash(group)
        T1 = T2 = T3 = T4 = D = R = 1
        X = pair(SK["D"], SK["D_hat"]) ** sk_s
        X = H.hashToZr(X)
        X = SK["g_hat"] ** X
        for name in Delta:
            tmp1 = SK[name]["T1"]
            tmp2 = SK[name]["T2"]
            tmp3 = SK[name]["T3"]
            tmp4 = SK[name]["T4"]

            T = Delta.copy()
            T.remove(name)
            for item in T:
                tmp1 *= SK[name]["Q"][item]["Q1"]
                tmp2 *= SK[name]["Q"][item]["Q2"]
                tmp3 *= SK[name]["Q"][item]["Q3"]
                tmp4 *= SK[name]["Q"][item]["Q4"]

            T1 *= tmp1
            T2 *= tmp2
            T3 *= tmp3
            T4 *= tmp4
            D *= SK[name]["D"]
            R *= (SK[name]["R"] / X)
            #R *= (SK[name]["R"])
        
        Y = pair(CT["C1"], D) * pair(CT["C2"], R) * pair(CT["E1"], T1) * pair(CT["E2"], T2) * pair(CT["E3"], T3) * pair(CT["E4"], T4)
        #return Y == (pk["Omega"] ** CT["mu"])
        return CT["C"] == Y
        
            

def main():
    #Get the eliptic curve with the bilinear mapping feature needed.
    #groupObj = PairingGroup('MNT224')
    groupObj = PairingGroup('BN254')
    Policy_Matrix = {'School':[1, 1, 0], "Dept":[0, -1, 1], "Deg":[0, 0, -1], "Pos":[0, -1, 0]}
    Policy = {'School':"NSYSU", "Dept":"CSE", "Deg":"MS", "Pos":"Teacher"}
    Keyword = {"School":"NSYSU", "Pos":"Teacher"}
    Delta = ["School", "Pos"]

    
    test_time = 100

    kpabks = MZNLHS17(groupObj)

    (msk, pk) = kpabks.setup()
    (pk_s, sk_s) = kpabks.s_keygen(pk)

    #CT = kpabks.encrypt(pk, Keyword)
    #SK = kpabks.keygen(msk, pk_s, Policy, Policy_Matrix)
    #tmp = kpabks.decrypt(pk, sk_s, SK, CT, Delta)
    #print(SK)

    start = time.time()
    for i in range(test_time):
        CT = kpabks.encrypt(pk, Keyword)
    end = time.time()
    print("Encryption: ", (end - start)/test_time)


    SK = kpabks.keygen(msk, pk_s, Policy, Policy_Matrix)

    start = time.time()
    for i in range(test_time):
        tmp = kpabks.decrypt(pk, sk_s, SK, CT, Delta)
    end = time.time()
    print("Decryption: ", (end - start)/test_time)
    print(tmp)
    fo = open("MZNLHS17_KeySize", "w")
    print(SK, file = fo)
    fo.close()
    fo = open("MZNLHS17_CiphertextSize", "w")
    print(CT, file = fo)
    fo.close()

if __name__ == '__main__':
    debug = True
    main()
