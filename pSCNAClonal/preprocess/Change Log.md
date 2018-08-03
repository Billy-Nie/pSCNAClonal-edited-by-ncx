# pSCNAClonal的修改日志qwq
<b>by 聂晨晞 2018年7月13日</b>
<br>修改目的:将原来phy-SCNAClonal中的向量操作修改为只支持一个的标量操作ヾ(◍°∇°◍)ﾉﾞ
<br>* converter.py中第180行,从原本的segPool修改为self._segPool
<br>* 132和133行,分别从原来的constants.STRIPENUML,constants.NOISESTRIPENUML修改为了constants.STRIPENUM和constants.NOISESTRIPENUM