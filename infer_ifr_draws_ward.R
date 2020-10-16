# adjusted Ward estimates

Sg=21.5 #55+ Rosenberg
w1=1008856 #55-64 pop
w2=689816 #65-74 pop
w3=845025 #75+ pop
s1=5.9 #ward 55-64
s2=3.2 #ward 65-74
s3=3.3 #ward 75+
# se2=(w1+w2)*Sg/(w1*(s1/s2)+w2)
# se1=(s1/s2)*se2

se1=(w1+w2+w3)*Sg/( w1+w2*s2/s1+w3*s3/s1)
se2=se1*s2/s1
se3=se1*s3/s1

# Scaling to IFRs from model

se1/se3*5.3
se1/se3*6.3

deaths_55to64 = 3556*(w1/2112562)
adj_45to64 = ((w1*se1 + 1103706*26.5)/2112562)/100

3556/(2112562*adj_45to64)*100 #45-64 IFR
1858/(1103706*26.5/100)*100 #45-54 IFR
deaths_55to64/(w1*se1/100)*100 #55-64 IFR
3963/(w2*se2/100)*100 #65-74 IFR
7731/(w3*se3/100)*100 #75+ IFR


