#!/bin/bash
for num in I1_sjureiXIII I2_sjureiXIII II III IV V VI1_sjureiVII VI2_sjureiVII VIII IX X XI XII1 XII2 XIV XV XVI MT; do
./mp_inference.py -n 150 HY3_sjurei/HY3_Hma/sjurei$num.txt HY3_sjurei/HY3_Lma/sjurei$num.txt -m contrast -o HY3_MAL_$num.out --plotFile HY3_MAL_$num |& tee HY3_MAL_$num.log
done

