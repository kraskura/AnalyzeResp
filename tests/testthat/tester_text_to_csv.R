# context("tester package here")

test_that(" This is used to test a function - make one for each function ", {

  d<-read.delim("./test_file_TESTER_MMR.txt", skip=19)

  expect_equal(d[4,15], 12.082) # to test function and make sure it cannot be easily "broken"

  expect_equal(d[16,5],8.0174)

})
