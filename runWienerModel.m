clear all; close all;

for fnum = 62:100
    fprintf('%i\n',fnum);
    wienerModel(fnum,'save');
end

!echo "laptop done" > msg.txt
!mail -s "SIM UPDATE" iobeid@vtext.com < msg.txt
!rm msg.txt
