'''
Yi-Fan Tseng, Chun-I Fan, Zi-Cheng Liu
 
| From: "Fast Keyword Search Over Encrypted Data with Short Ciphertext in Clouds"
| Published in: 
| Available from: 
| Notes: Security Assumption: DBDH on asymmetric paring group. 

* type:           SE from KPABE
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
class TFL19(ABEnc):
    
    def __init__(self, groupObj):
        ABEnc.__init__(self)
        global util, group
        util = SecretUtil(groupObj, debug)        
        group = groupObj

    def setup(self):
        g1, g2 = group.random(G1), group.random(G2)
        alpha, beta, phi = group.random(), group.random(), group.random()        
        h = g1 ** phi
        h_hat = g2 ** phi
        U = pair(g1, g2) ** (alpha * (beta - 1))
        V = pair(g1, g2) ** (alpha * beta)
               
        pk = {'g1':g1, 'g2':g2, 'U':U, 'V':V, 'h':h}
        msk = pk.copy()
        msk["g_hat^alpha"] = g2 ** alpha
        msk["h_hat"] = h_hat
        return (msk, pk)

    def encrypt(self, pk, M, Attribute):
        H = Hash(group)
        k = group.random()

        CT = dict()
        CT["C1"] = M * pk["V"]** k
        CT["C2"] = pk["U"] ** k
        CT["C3"] = pk["g1"] ** k

        for name, value in Attribute.items():
            r = H.hashToZr(name + value)
            CT[name] = (pk["h"] * pk["g1"] ** r)**k
        
        return CT
    
    def keygen(self, pk, msk, Pol, Pol_M):
        H = Hash(group)
        SK = dict()
        s, y2, y3 = 1, group.random(), group.random()

        sigma = dict()
        for name, value in Pol.items():
            h = H.hashToZr(name + value)
            sigma[name] = msk["h_hat"] * msk["g2"]**h
        
        for name, row in Pol_M.items():
            lamb = row[0] + row[1]*y2 + row[2]*y3
            r = group.random()
            
            Q = dict()
            D0 = (msk["g_hat^alpha"]**lamb) * (sigma[name] ** r)
            D1 = msk["g2"]**r
            tmp = Pol.copy()
            del tmp[name]
            for j in tmp:
                Q[j] = sigma[j] ** r
            SK[name] = {"D0":D0, "D1":D1, "Q":Q}

        return SK
    
    def decrypt(self, pk, SK, CT, Delta):
        
        L = numerator = denominator = 1
        for name in Delta:
            L *= CT[name]
        for name in Delta:
            tmp = SK[name]["D0"]
            T = Delta.copy()
            T.remove(name)
            for item in T:
                tmp *= SK[name]["Q"][item]
            numerator *= tmp
            denominator *= SK[name]["D1"]
        Z = (pair(CT["C3"], numerator)) / (pair(L, denominator))
        m = CT["C1"] / (CT["C2"] * Z)
        return m
        
            

def main():
    #Get the eliptic curve with the bilinear mapping feature needed.
    #groupObj = PairingGroup('MNT224')
    groupObj = PairingGroup('BN254')
    Policy_Matrix = {'School':[1, 1, 0], "Dept":[0, -1, 1], "Deg":[0, 0, -1], "Pos":[0, -1, 0]}
    Policy = {'School':"NSYSU", "Dept":"CSE", "Deg":"MS", "Pos":"Teacher"}
    Attribute = {"School":"NSYSU", "Pos":"Teacher"}
    Delta = ["School", "Pos"]

    test_time = 100

    kpabe = TFL19(groupObj)

    (msk, pk) = kpabe.setup()

    m = group.random(GT)

    start = time.time()
    for i in range(test_time):
        CT = kpabe.encrypt(pk, m, Attribute)
    end = time.time()
    print("Encryption: ", (end - start)/test_time)


    SK = kpabe.keygen(pk, msk, Policy, Policy_Matrix)

    start = time.time()
    for i in range(test_time):
        msg = kpabe.decrypt(pk, SK, CT, Delta)
    end = time.time()
    print("Decryption: ", (end - start)/test_time)

    #print(m == msg)
    fo = open("TFL19_KeySize", "w")
    print(SK, file = fo)
    fo.close()
    fo = open("TFL19_CiphertextSize", "w")
    print(CT, file = fo)
    fo.close()

if __name__ == '__main__':
    debug = True
    main()
