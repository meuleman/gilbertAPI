To start the SSL-backed server:

```
$ PORT=9877 node server.js
```

To test a request from another host:

```
$ args=(-d '--header=reqState:chr11:44556463-44572536' '--header=purpose:getRegions' '--header=factors:[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]' '--header=detailLevel:1')
$ wget -qO- "${args[@]}" https://tools.altiusinstitute.org:9877 | jq .
...
```

Replace 9877 with whichever port is publicly open or otherwise available between server and client.
