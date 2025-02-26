import sys
import os
import time

if __name__ == "__main__":
    usage = "python %s <old_run_case> <new_run_case>\n" % sys.argv[0]
    usage += "e.g. python %s " % sys.argv[0]
    usage += "zjet2d-rlc50-sig10-delta0p1-vinjvd-rotramp "
    usage += "zjet2d-rlc50-sig10-delta0p1-vj0p5-rotramp"
    if len(sys.argv) != 3:
        print(usage)
        sys.exit(0)
    old_test_case_dir, new_test_case_dir = sys.argv[1:]

    old_test_case = os.path.basename(os.path.normpath(old_test_case_dir))
    new_test_case = os.path.basename(os.path.normpath(new_test_case_dir))

    dirfiles = os.listdir(old_test_case_dir)
    copy_these = []
    copy_these_to = []
    for dirfile in dirfiles:
        if old_test_case in dirfile and not dirfile.startswith("."):
            copiedfile = os.path.join(old_test_case_dir, dirfile) 
            copy_these.append(copiedfile)
            newfile = os.path.join(new_test_case_dir,
                dirfile.replace(old_test_case, new_test_case))
            copy_these_to.append(newfile)

    if os.path.isdir(new_test_case_dir):
        print("About to make the following copies")
    else:
        print("About to create directory: %s" % new_test_case_dir)
        time.sleep(2)
        print("and also make the following copies...")
    time.sleep(2)

    for cp, cpto in zip(copy_these, copy_these_to):
        print("\t%s --> %s" % (cp, cpto))
    time.sleep(2)
    print("Proceed? Enter y/n... ")

    try:
        result = raw_input()
    except:
        result = input()
    should_proceed = result.lower().startswith("y")

    if should_proceed:
        if not os.path.isdir(new_test_case_dir):
            try:
                os.mkdir(new_test_case_dir)
            except:
                print("Cannot make directory %s")
                print("A parent directory in the path does not exist.")
                print("Please create parent directories and try again.")
                sys.exit(0)
            
        print("Copying files...")
        for cp, cpto in zip(copy_these, copy_these_to):
            # Don't do this; I need to replace occurrences of old_test_case
            # with new_test_case in each file
            # shutil.copyfile(cp, cpto)
            srcfile = open(cp, 'r')
            lines = srcfile.readlines()
            srcfile.close()
            dstfile = open(cpto, 'w')
            for line in lines:
                dstfile.write(line.replace(old_test_case, new_test_case))
            dstfile.close()
        print("Done.")
    else:
        print("Operation aborted.")
