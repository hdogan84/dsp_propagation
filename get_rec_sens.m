function [rec_out]=get_rec_sens(fin,rec);
% Function which outputs hydrophone sensitivity (rec_out) for 
% input frequencies (fin in Hz) for input hydrophone serial 
% number (rec: use last two digits only). Sensitivities also 
% available in equipment folder
if rec==34;
    rec_sens=[1500	-209.2
2000	-208.3
2500	-209.9
3000	-210.4
3500	-209
4000	-209.4
4500	-208.9
5000	-209.2
5500	-208.7
6000	-208.7
6500	-208.4
7000	-208.7
7500	-208.8
8000	-208.6
8500	-209.3
9000	-209.1
9500	-209.3
10000	-209.7
11000	-210.1
12000	-209.7
13000	-209.5
14000	-209.7
15000	-209.9
16000	-210.1
17000	-210.5
18000	-210.4
19000	-210.8
20000	-211.4
21000	-211.1
22000	-210.8
23000	-210.7
24000	-210.3
25000	-210.6
26000	-210.5
27000	-210.8
28000	-210.4
29000	-210.5
30000	-210.3
31000	-210.5
32000	-210.8
33000	-211
34000	-211.4
35000	-211.5
36000	-211.4
37000	-211.3
38000	-211.4
39000	-211.3
40000	-211.2
41000	-211.4
42000	-211.6
43000	-211.8
44000	-211.8
45000	-211.9
46000	-211.9
47000	-211.9
48000	-211.8
49000	-211.8
50000	-211.9
51000	-212
52000	-212.1
53000	-212.3
54000	-212.5
55000	-212.6
56000	-212.8
57000	-213
58000	-213.3
59000	-213.4
60000	-213.6
61000	-213.8
62000	-213.9
63000	-214.3
64000	-214.8
65000	-215.2
66000	-215.6
67000	-215.9
68000	-216.1
69000	-215.7
70000	-215.4
71000	-215.3
72000	-215.1
73000	-214.9
74000	-214.8
75000	-214.6
76000	-214.4
77000	-214.4
78000	-214.5
79000	-214.5
80000	-214.6
81000	-214.9
82000	-214.9
83000	-214.8
84000	-214.8
85000	-214.8
86000	-215
87000	-215
88000	-215
89000	-215
90000	-215
91000	-215
92000	-215.1
93000	-215.1
94000	-215.2
95000	-215.3
96000	-215.4
97000	-215.6
98000	-215.6
99000	-215.6
100000	-215.8
101000	-216
102000	-216.3
103000	-216.4
104000	-216.6
105000	-216.6
106000	-216.6
107000	-216.5
108000	-216.5
109000	-216.5
110000	-216.4
111000	-216.2
112000	-215.9
113000	-215.8
114000	-215.6
115000	-215.4
116000	-215.3
117000	-215
118000	-214.7
119000	-214.5
120000	-214.3
121000	-214.1
122000	-213.8
123000	-213.5
124000	-213.3
125000	-213.1
130000	-212.4
135000	-212.7
140000	-213.3
145000	-211.4
150000	-208.8
155000	-206.4
160000	-204
165000	-202.4
170000	-202.4
175000	-204.3
180000	-206.7
185000	-209.5
190000	-211.9
195000	-213.7
200000	-215.1
205000	-216
210000	-217.6];
elseif rec==38;
    rec_sens=[2000	-210.1
2500	-209.2
3000	-210.5
3500	-208.8
4000	-208.8
4500	-208.9
5000	-209.2
5500	-208.7
6000	-208.6
6500	-208.1
7000	-208.3
7500	-208.4
8000	-208.4
8500	-208.9
9000	-209.1
9500	-209.2
10000	-209.6
11000	-210
12000	-209.9
13000	-209.8
14000	-210
15000	-210.4
16000	-210.4
17000	-210.4
18000	-210.1
19000	-210.6
20000	-210.3
21000	-210.1
22000	-209.5
23000	-208.9
24000	-209.1
25000	-209.4
26000	-209.3
27000	-209.6
28000	-209.1
29000	-209
30000	-209.1
31000	-209.3
32000	-209.6
33000	-210.1
34000	-210.6
35000	-210.8
36000	-210.9
37000	-210.9
38000	-210.9
39000	-211.1
40000	-211.5
41000	-211.9
42000	-212.1
43000	-212.1
44000	-212.2
45000	-212.4
46000	-212.6
47000	-212.7
48000	-212.5
49000	-212.5
50000	-212.5
51000	-212.7
52000	-212.8
53000	-212.9
54000	-212.9
55000	-212.8
56000	-213
57000	-213.2
58000	-213.3
59000	-213.2
60000	-213.3
61000	-213.2
62000	-213
63000	-212.9
64000	-212.9
65000	-212.8
66000	-212.7
67000	-212.8
68000	-213.2
69000	-213.2
70000	-213.3
71000	-213.6
72000	-213.8
73000	-213.9
74000	-214.1
75000	-214
76000	-213.9
77000	-213.8
78000	-213.9
79000	-213.7
80000	-214.1
81000	-214.4
82000	-214.3
83000	-214.3
84000	-214.2
85000	-214
86000	-214.1
87000	-214
88000	-213.8
89000	-213.7
90000	-213.6
91000	-213.5
92000	-213.6
93000	-213.5
94000	-213.4
95000	-213.5
96000	-213.4
97000	-213.5
98000	-213.5
99000	-213.4
100000	-213.2
101000	-213.2
102000	-213.1
103000	-213
104000	-212.9
105000	-212.7
106000	-212.5
107000	-212.4
108000	-212.3
109000	-212.2
110000	-212.2
111000	-211.9
112000	-211.7
113000	-211.6
114000	-211.6
115000	-211.4
116000	-211.3
117000	-211.2
118000	-211
119000	-210.9
120000	-210.7
121000	-210.7
122000	-210.5
123000	-210.5
124000	-210.1
125000	-209.9
130000	-209.5
135000	-210.1
140000	-208.8
145000	-207.2
150000	-205.2
155000	-203.5
160000	-202
165000	-200.9
170000	-201.2
175000	-203.4
180000	-205.6
185000	-207.2
190000	-209
195000	-210.6
200000	-212
205000	-212.9
210000	-214.3];
elseif rec==39;
    rec_sens=[2000	-209.9
2500	-209.8
3000	-210.7
3500	-209.1
4000	-208.8
4500	-208.9
5000	-208.6
5500	-208.9
6000	-208.7
6500	-208.2
7000	-208.7
7500	-208.8
8000	-208.8
8500	-209.6
9000	-209.6
9500	-209.6
10000	-210
11000	-210.6
12000	-210.5
13000	-210
14000	-210.2
15000	-210.4
16000	-210.5
17000	-210.6
18000	-210.1
19000	-210.6
20000	-210.8
21000	-210.7
22000	-210
23000	-209.4
24000	-209.5
25000	-209.8
26000	-209.9
27000	-210.3
28000	-210.1
29000	-210.1
30000	-210
31000	-210
32000	-210.3
33000	-210.8
34000	-211.3
35000	-211.4
36000	-211.3
37000	-211.4
38000	-211.5
39000	-211.7
40000	-212
41000	-212.4
42000	-212.7
43000	-212.7
44000	-212.7
45000	-212.9
46000	-213.1
47000	-213.1
48000	-213
49000	-212.9
50000	-213
51000	-213.1
52000	-213.2
53000	-213.3
54000	-213.2
55000	-213.1
56000	-213.1
57000	-213.1
58000	-213.1
59000	-212.9
60000	-212.9
61000	-212.7
62000	-212.5
63000	-212.4
64000	-212.4
65000	-212.3
66000	-212.2
67000	-212.2
68000	-212.2
69000	-212.1
70000	-212.1
71000	-212.5
72000	-212.6
73000	-212.7
74000	-212.7
75000	-212.6
76000	-212.5
77000	-212.5
78000	-212.5
79000	-212.5
80000	-212.8
81000	-213
82000	-213
83000	-212.9
84000	-212.8
85000	-212.9
86000	-213
87000	-212.9
88000	-212.9
89000	-212.8
90000	-212.8
91000	-212.8
92000	-212.9
93000	-213
94000	-213.1
95000	-213.3
96000	-213.5
97000	-213.6
98000	-213.7
99000	-213.7
100000	-213.6
101000	-213.5
102000	-213.2
103000	-213.1
104000	-213
105000	-212.8
106000	-212.8
107000	-212.6
108000	-212.3
109000	-212.4
110000	-212.4
111000	-212.2
112000	-212
113000	-211.9
114000	-211.8
115000	-211.8
116000	-211.7
117000	-211.5
118000	-211.3
119000	-211.2
120000	-211.1
121000	-211
122000	-210.9
123000	-210.7
124000	-210.5
125000	-210.5
130000	-209.7
135000	-209.7
140000	-209.7
145000	-208.8
150000	-207.4
155000	-206.2
160000	-204.6
165000	-203.9
170000	-203.5
175000	-204.2
180000	-205.8
185000	-207.7
190000	-210.2
195000	-212.9
200000	-214.6
205000	-215.9
210000	-217.9];
elseif rec==36;
    rec_sens=[2000	-207.5
2500	-210.1
3000	-210.5
3500	-209.7
4000	-209.3
4500	-209.6
5000	-208.8
5500	-209.2
6000	-208.9
6500	-208.4
7000	-208.9
7500	-208.8
8000	-208.6
8500	-209.4
9000	-209.1
9500	-209.3
10000	-209.7
11000	-210.3
12000	-209.8
13000	-209.5
14000	-209.6
15000	-209.8
16000	-210.1
17000	-210.6
18000	-210.4
19000	-210.9
20000	-211
21000	-211.1
22000	-210.5
23000	-209.9
24000	-209.6
25000	-209.9
26000	-209.8
27000	-210.2
28000	-209.9
29000	-209.9
30000	-209.9
31000	-210
32000	-210.2
33000	-210.6
34000	-211.2
35000	-211.3
36000	-211.2
37000	-211
38000	-211
39000	-211.1
40000	-211.1
41000	-211.3
42000	-211.5
43000	-211.7
44000	-211.8
45000	-211.9
46000	-212
47000	-212
48000	-211.9
49000	-211.8
50000	-211.8
51000	-211.8
52000	-211.9
53000	-212
54000	-212.1
55000	-212
56000	-212.1
57000	-212.3
58000	-212.5
59000	-212.6
60000	-212.7
61000	-212.8
62000	-212.8
63000	-212.8
64000	-212.8
65000	-212.8
66000	-212.6
67000	-212.6
68000	-212.6
69000	-212.4
70000	-212.3
71000	-212.6
72000	-212.7
73000	-212.7
74000	-212.8
75000	-212.8
76000	-212.7
77000	-212.7
78000	-212.9
79000	-213
80000	-213.3
81000	-213.7
82000	-213.6
83000	-213.5
84000	-213.5
85000	-213.5
86000	-213.7
87000	-213.8
88000	-213.7
89000	-213.7
90000	-213.5
91000	-213.5
92000	-213.5
93000	-213.4
94000	-213.2
95000	-213.2
96000	-212.9
97000	-212.7
98000	-212.8
99000	-212.8
100000	-212.8
101000	-212.9
102000	-212.9
103000	-213.1
104000	-213
105000	-212.7
106000	-212.7
107000	-212.6
108000	-212.6
109000	-212.6
110000	-212.6
111000	-212.4
112000	-212.2
113000	-212.1
114000	-212
115000	-211.9
116000	-211.8
117000	-211.5
118000	-211.3
119000	-211.1
120000	-210.9
121000	-210.7
122000	-210.5
123000	-210.3
124000	-210
125000	-210
130000	-209.6
135000	-210.2
140000	-209.7
145000	-208.5
150000	-207.2
155000	-205.8
160000	-204.8
165000	-203.9
170000	-204
175000	-205.4
180000	-206.4
185000	-207.9
190000	-209.5
195000	-211.5
200000	-213.1
205000	-214.1
210000	-215.6];
elseif rec==37;
    rec_sens=[
2000	-207.6
2500	-210.5
3000	-209.5
3500	-209.2
4000	-208.5
4500	-208.7
5000	-208.5
5500	-208.3
6000	-208.4
6500	-207.7
7000	-208.1
7500	-208.2
8000	-208.1
8500	-208.8
9000	-208.8
9500	-208.7
10000	-209.1
11000	-209.5
12000	-209.4
13000	-209.2
14000	-209.1
15000	-209.2
16000	-209.3
17000	-209.8
18000	-209.9
19000	-210
20000	-210.4
21000	-210.4
22000	-209.9
23000	-209.5
24000	-209.3
25000	-209.4
26000	-209.2
27000	-209.8
28000	-209.6
29000	-209.7
30000	-209.7
31000	-209.8
32000	-209.9
33000	-210.1
34000	-210.8
35000	-210.9
36000	-210.6
37000	-210.3
38000	-210.4
39000	-210.7
40000	-210.8
41000	-210.9
42000	-211.2
43000	-211.5
44000	-211.6
45000	-211.8
46000	-211.8
47000	-211.8
48000	-211.6
49000	-211.5
50000	-211.5
51000	-211.5
52000	-211.4
53000	-211.5
54000	-211.6
55000	-211.5
56000	-211.7
57000	-211.8
58000	-212
59000	-212.1
60000	-212.3
61000	-212.4
62000	-212.3
63000	-212.4
64000	-212.5
65000	-212.4
66000	-212.2
67000	-212.2
68000	-212.1
69000	-212
70000	-211.9
71000	-212
72000	-212.1
73000	-212.1
74000	-212.1
75000	-212
76000	-211.9
77000	-212
78000	-212.1
79000	-212.1
80000	-212.5
81000	-212.8
82000	-212.8
83000	-212.7
84000	-212.6
85000	-212.6
86000	-212.6
87000	-212.6
88000	-212.5
89000	-212.4
90000	-212.3
91000	-212.3
92000	-212.3
93000	-212.3
94000	-212.2
95000	-212.3
96000	-212.5
97000	-212.4
98000	-212.4
99000	-212.3
100000	-212.2
101000	-212.2
102000	-212.1
103000	-212.2
104000	-212.1
105000	-212
106000	-212
107000	-211.9
108000	-211.8
109000	-211.8
110000	-211.8
111000	-211.6
112000	-211.4
113000	-211.1
114000	-211.1
115000	-211.1
116000	-211.1
117000	-210.9
118000	-210.7
119000	-210.6
120000	-210.5
121000	-210.4
122000	-210.2
123000	-209.9
124000	-209.8
125000	-209.6
130000	-209.1
135000	-209.5
140000	-209.8
145000	-208.4
150000	-207.6
155000	-206.2
160000	-205
165000	-204.1
170000	-203.9
175000	-205
180000	-205.9
185000	-206.8
190000	-208.3
195000	-210.5
200000	-212.2
205000	-213.1
210000	-214.7];
else xxx
end

% interpolate receiver sensitivity
%figure;plot(rec_sens(:,1),rec_sens(:,2),'.')
finterp=[2e3:500:210e3];
rec=interp1(rec_sens(:,1),rec_sens(:,2),finterp,'spline');
%hold on;plot(finterp,rec,'r');
rec_out=zeros(size(fin));
 for n=1:length(fin);
     f=fin(n);
     find(finterp==f);
     rec_out(n)=rec(find(finterp==f));
 end
 %plot(fin,rec_out,'go');