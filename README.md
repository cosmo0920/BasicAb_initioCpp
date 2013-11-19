basic Ab initio Cpp
===

HeH+の第一原理計算を行うアプリケーション。
第一原理計算の再初歩。

2電子分子ならパラメーターを変えれば計算できる（と思います）。

3電子以上は計算の方法を変える必要があります。が、そこまでは書けていません。

このプログラムは hehp.f [code](http://www.chemie.unibas.ch/~meuwly/pdfs/hehp.f) をC++に翻訳したものです。

計算のアルゴリズムの勉強用にかなり適当に組んでいるのでC++らしくないところが多々あります。

## COPYING

### license

LGPL 2.1.

ただしclapack.hは除く。

see: [COPYING.LIB](COPYING.LIB)