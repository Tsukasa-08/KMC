# (1)コンパイラ
CC  = g++
# (2)コンパイルオプション
CFLAGS    =
# (3)実行ファイル名
TARGET  = test_kmc_test
# (4)コンパイル対象のソースコード
SRCS    = ./kmc_test.cpp
# (5)オブジェクトファイル名
OBJS    = $(SRCS:.cpp=.o)
 
# (6)インクルードファイルのあるディレクトリパス
# 必要なものはHOMEディレクトリ以下のincludeディレクトリに入れとく
INCDIR  = -I../inc
INCDIR  += -I$(HOME)/include/eigen-3.3.9
INCDIR  += -I$(HOME)/include
 
# (7)ライブラリファイルのあるディレクトリパス
LIBDIR  = 
 
# (8)追加するライブラリファイル
LIBS    = 

# (9)ターゲットファイル生成
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(LIBDIR) $(LIBS)
	
# (10)オブジェクトファイル生成
$(OBJS): $(SRCS)
	$(CC) $(CFLAGS) $(INCDIR) -c $(SRCS)

# (11)"make all"で make cleanとmakeを同時に実施。
all: clean $(OBJS) $(TARGET)
# (12).oファイル、実行ファイル、.dファイルを削除
clean:
	-rm -f $(OBJS) $(TARGET) *.d

#tsukasaが追加、作成したプログラムを実行する
run:
	./do_kmc.sh
