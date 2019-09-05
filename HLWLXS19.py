'''
Jialu Hao, Jian Liu, Huimei Wang, Lingshuang Liu, Ming Xian, Xuemin Shen
 
| From: "Efficient Attirubte-Based Access Control with Authorized Search in Cloud Storage"
| Published in: 
| Available from: 
| Notes: Security Assumption:

* type:           SE supporting monotonic search query
* setting:        Pairing

:Authors:    Yi-Fan Tseng
:Date:            09/03/2019
'''
from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from charm.toolbox.secretutil import SecretUtil
from charm.toolbox.ABEnc import ABEnc
from charm.toolbox.hash_module import Waters, Hash
import time

debug = False
class HLWLXS19(ABEnc):
    
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
        w = g ** x3
        w_hat = g_hat ** x3
        #u, h, w= group.random(G1), group.random(G1), group.random(G1)
        #u_hat, h_hat, w_hat = group.random(G2), group.random(G2), group.random(G2)
        alpha, tau1, tau2, tau3, tau4 = group.random(), group.random(), group.random(), group.random(), group.random()
        g1 = g**tau1
        g2 = g**tau2
        g3 = g**tau3
        g4 = g**tau4
        pk = {"g":g, "u":u, "h":h, "w":w, "g1":g1, "g2":g2, "g3":g3, "g4":g4, "Omega":pair(g, g_hat)**alpha, "omega":pair(g, g_hat)}
        msk = pk.copy()
        msk["alpha"] = alpha
        msk["g_hat"] = g_hat
        msk["u_hat"] = u_hat
        msk["h_hat"] = h_hat
        msk["w_hat"] = w_hat
        msk["tau1"] = tau1
        msk["tau2"] = tau2
        msk["tau3"] = tau3
        msk["tau4"] = tau4

        return (msk, pk)

    def encrypt(self, pk, Keyword, M):
        H = Hash(group)
        s = group.random()

        CT = dict()
        E_tilde = pk["Omega"]**s
        E_tilde *= M
        E = pk["g"]**s
        CT["E_tilde"] = E_tilde
        CT["E"] = E

        for name, value in Keyword.items():
            z, s1, s2 = group.random(), group.random(), group.random()
            htemp = H.hashToZr(name + value)
            
            E0 = (pk["w"] ** -s) * ( ((pk["u"] ** htemp) * pk["h"]) ** z )
            E1 = pk["g1"] ** (z - s1)
            E2 = pk["g2"] ** s1
            E3 = pk["g3"] ** (z - s2)
            E4 = pk["g4"] ** s2
            CT[name] = {"E0":E0, "E1":E1, "E2":E2, "E3":E3, "E4":E4}
        #CT["s"] = s
        return CT
    
    def keygen(self, msk, Pol):
        H = Hash(group)
        SK = dict()
        
        tau1 = msk["tau1"]
        tau2 = msk["tau2"]
        tau3 = msk["tau3"]
        tau4 = msk["tau4"]

        '''
        The access structure is fixed for "NSYSU" AND ( ("CSE" AND "Mas") OR "Teacher" ). This implementation is only for efficiency test.
        Label: Root -> 0, NSYSU -> 1, OR GATE -> 2, AND GATE -> 3, Teacher -> 4, CSE -> 5, Mas -> 6
        '''
        leaf = ["NSYSU", "Teacher", "CSE", "Mas"]
        #q0(x) = alpha + q0*x
        q0 = group.random()
        #q1 = q0(1)
        q1 = msk["alpha"] + q0 
        #q2 = q0(2)
        q4 = q2 = msk["alpha"] + q0*2
        #q3(x) = q2(3) + q3*x = q2 + q3*x
        q3 = group.random()
        #q4 = q2(4) = q0(2)
        #q5 = q3(5)
        q5 = q2 + q3 * 5
        #q6 = q3(6)
        q6 = q2 + q3 * 6
        q = dict()
        q["NSYSU"] = q1
        q["CSE"] = q5
        q["Mas"] = q6
        q["Teacher"] = q4
                
        
        for name, value in Pol.items():
            t1, t2 = group.random(), group.random()

            t = tau1*tau2*t1 + tau3*tau4*t2
            h = H.hashToZr(name + value)
            Y = (msk["u_hat"]**h) * msk["h_hat"]

            D = ( msk["g_hat"] ** q[value] ) * ( msk["w_hat"] **  t)
            D0 = msk["g_hat"]** t
            D1 = Y ** (-tau2*t1)
            D2 = Y ** (-tau1*t1)
            D3 = Y ** (-tau4*t2)
            D4 = Y ** (-tau3*t2)

            SK[name] = {"D":D, "D0":D0, "D1":D1, "D2":D2, "D3":D3, "D4":D4}
        #SK["q"] = q
        #SK["alpha"] = msk["alpha"]

        return SK

    def decrypt(self, pk, SK, CT, Delta):
        H = Hash(group)
        
        P_School = ( pair(CT["E"], SK["School"]["D"]) * pair(CT["School"]["E0"], SK["School"]["D0"]) * pair(CT["School"]["E1"], SK["School"]["D1"])
                 * pair(CT["School"]["E2"], SK["School"]["D2"]) * pair(CT["School"]["E3"], SK["School"]["D3"]) * pair(CT["School"]["E4"], SK["School"]["D4"]) )
        P_Pos = ( pair(CT["E"], SK["Pos"]["D"]) * pair(CT["Pos"]["E0"], SK["Pos"]["D0"]) * pair(CT["Pos"]["E1"], SK["Pos"]["D1"])
                 * pair(CT["Pos"]["E2"], SK["Pos"]["D2"]) * pair(CT["Pos"]["E3"], SK["Pos"]["D3"]) * pair(CT["Pos"]["E4"], SK["Pos"]["D4"]) )
        
        #s = CT["s"]
        #q = SK["q"]
        #alpha = SK["alpha"]
        #print(2* q["NSYSU"] - q["Teacher"] == alpha)
        #print( pk["omega"]**(q["NSYSU"]*s) == P_School)

        K = (P_School ** 2 ) * ( P_Pos ** (-1) )
        
        M = CT["E_tilde"] / K
        
        return M
        
            

def main():
    #Get the eliptic curve with the bilinear mapping feature needed.
    groupObj = PairingGroup('MNT224')
    #groupObj = PairingGroup('BN254')
    Policy_Matrix = {'School':[1, 1, 0], "Dept":[0, -1, 1], "Deg":[0, 0, -1], "Pos":[0, -1, 0]}
    Policy = {'School':"NSYSU", "Dept":"CSE", "Deg":"Mas", "Pos":"Teacher"}
    Keyword = {"School":"NSYSU", "Pos":"Teacher"}
    Delta = ["School", "Pos"]

    
    test_time = 100

    kpabe = HLWLXS19(groupObj)

    (msk, pk) = kpabe.setup()

    M = group.random(GT)
    #CT = kpabe.encrypt(pk, Keyword, M)
    #SK = kpabe.keygen(msk, Policy)
    #print(2*SK["NSYSU"]-SK["Teacher"] == msk["alpha"])
    #m = kpabe.decrypt(pk, SK, CT, Delta)
    #print(M == m)
    
    
    start = time.time()
    for i in range(test_time):
        CT = kpabe.encrypt(pk, Keyword, M)
    end = time.time()
    print("Encryption: ", (end - start)/test_time)


    SK = kpabe.keygen(msk, Policy)
    
    start = time.time()
    for i in range(test_time):
        m = kpabe.decrypt(pk, SK, CT, Delta)
    end = time.time()
    print("Decryption: ", (end - start)/test_time)
    
    print(m ==M)
    fo = open("HLWLXS19_KeySize", "w")
    print(SK, file = fo)
    fo.close()
    fo = open("HLWLXS19_CiphertextSize", "w")
    print(CT, file = fo)
    fo.close()
    
    
    

if __name__ == '__main__':
    debug = True
    main()
