class create_starting_sample:
    def get_sample_manual(self,epi,depth,strike,dip,rake,M0,savepath):
        with open(savepath, 'w') as save_file:
            save_file.write("%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n\r" % (epi,depth,strike,dip,rake,M0))
        save_file.close()
        print("Txt file with new sample is saved as %s" % savepath)
        return savepath