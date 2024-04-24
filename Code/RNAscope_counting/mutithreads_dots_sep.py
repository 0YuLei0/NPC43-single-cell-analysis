from ij import IJ, ImagePlus 
from ij.process import FloatProcessor 
import threading

def rm(x,n):

	for i in range(x,x+n):
		
        #remove the red dot
		if normredpixels[i] > normgreenpixels[i]:
			newredpixels[i]=255
	
        #remove the green dot
		if normgreenpixels[i] > normredpixels[i]: #-12 					
			newgreenpixels[i]=255
	
        #remove background for red dot
		if (normredpixels[i]+normgreenpixels[i]) > 320: 					
			newgreenpixels[i]=255
			#newredpixels[i]=255

		#remove background for green dot
		if (normredpixels[i]+normgreenpixels[i]) > 330: 					
			#newgreenpixels[i]=255
			newredpixels[i]=255

        #remove cell without removing overlapped dots for green dots only
		if normgreenpixels[i]-normredpixels[i] < 5: #and normgreenpixels[i]+normredpixels[i] > 200: 									
			newredpixels[i]=255
	
        #remove cell without removing overlapped dots for red dots only
		if abs(normgreenpixels[i]-normredpixels[i]) < 10: # and normgreenpixels[i]+normredpixels[i] > 250:					
			newgreenpixels[i]=255


imp = IJ.getImage()  
cp = imp.getProcessor()# as a copy  
red = cp.toFloat(0, None)  
green = cp.toFloat(1, None)
newred = red.duplicate()
newgreen = green.duplicate()
redpixels = red.getPixels()
greenpixels = green.getPixels()
redmean= sum(redpixels)/len(redpixels)#-77
greenmean= sum(greenpixels)/len(greenpixels)#-77
redratio=161/redmean
greenratio=160/greenmean
normredpixels=[i*redratio for i in redpixels]
normgreenpixels=[i*greenratio for i in greenpixels]
newredpixels = newred.getPixels() 
newgreenpixels = newgreen.getPixels() 
t=4
n=len(redpixels)//t

threads = list()	
for x in range(0,n*t,n):
	thread = threading.Thread(target = rm(x,n))
	thread.start()
	threads.append(thread)	
for thread in threads:
	thread.join()

#print("newredpixels",newredpixels)
#print("newgreenpixels",newgreenpixels)
ip2=FloatProcessor(cp.width, cp.height, newredpixels, None)
ip2=ip2.convertToByte(True) #convert to 8 bits
imp2=ImagePlus("greendots", ip2)  

ip3=FloatProcessor(cp.width, cp.height, newgreenpixels, None) 
ip3=ip3.convertToByte(True) #convert to 8 bits
imp3=ImagePlus("reddots", ip3)

imp2.show()  
imp3.show()  



