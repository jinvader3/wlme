chirp.npy:
	python3 ../chirp_zadoff_chu.py --nzc 353 --q 1 --u 7 --out chirp.npy

rx0.npy: testchamber0.wlmadata chirp.npy
	python3 ../meshcage2.py --tx chirp.npy --wlma testchamber0.wlmadata --wv 1490 --sps 256000 --time 1 --cycles 1000 --osps 256000 --outpath rx0.npy

image0.data: chirp.npy rx0.npy
	python3 ../imager2.py --wv 1490 --sps 256000 --tx chirp.npy --wlma testchamber0.wlmadata --rx rx0.npy --image image0.png --cx 0 --cy 0 --cz 0 --w 120 --r 200

rx1.npy: testchamber1.wlmadata chirp.npy
	python3 ../meshcage2.py --tx chirp.npy --wlma testchamber1.wlmadata --wv 1490 --sps 256000 --time 1 --cycles 1000 --osps 256000 --outpath rx1.npy

image1.data: chirp.npy rx1.npy
	python3 ../imager2.py --wv 1490 --sps 256000 --tx chirp.npy --wlma testchamber0.wlmadata --rx rx1.npy --image image1.png --cx 0 --cy 0 --cz 0 --w 120 --r 200

rx2.npy: testchamber2.wlmadata chirp.npy
	python3 ../meshcage2.py --tx chirp.npy --wlma testchamber2.wlmadata --wv 1490 --sps 256000 --time 1 --cycles 1000 --osps 256000 --outpath rx2.npy

image2.data: chirp.npy rx2.npy
	python3 ../imager2.py --wv 1490 --sps 256000 --tx chirp.npy --wlma testchamber2.wlmadata --rx rx2.npy --image image2.png --cx 0 --cy 0 --cz 0 --w 120 --r 200

all: image0.data image1.data image2.data
