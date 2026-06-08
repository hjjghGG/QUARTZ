@echo off
tshark -r equinix-nyc.dirA.20190117-130200.UTC.anon.pcap -T fields -e frame.time_relative -e ip.src -e ip.dst -e tcp.srcport -e tcp.dstport -e udp.srcport -e udp.dstport -E separator=" " > CAIDA.txt
echo complete tuesday.csv
pause