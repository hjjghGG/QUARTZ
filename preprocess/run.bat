@echo off
tshark -r 202412251400.pcap -T fields -e frame.time_relative -e ip.src -e ip.dst -e tcp.srcport -e tcp.dstport -e udp.srcport -e udp.dstport -E separator="," > MAWI_.csv
echo complete tuesday.csv
pause